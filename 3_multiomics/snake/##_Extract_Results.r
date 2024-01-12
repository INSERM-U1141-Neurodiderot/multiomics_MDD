library (dplyr)
library (MOFA2)

num_factors = 10 

# set in c('Female_MDD', 'Male_MDD')
# o in c('DNAm', 'mRNA', 'miRNA', ''

set = 'Female_MDD'
o = 'mRNA'

omics = readRDS (paste('01_Datasets/', set, '_omic_', o, '.RDS', sep = ''))

################
###   MCIA   ###
################
MCIA_dry = readRDS (paste('03_Results/', set, '/01_MCIA/', o, '/01_factorizations_mcia_dry.RDS', sep = ''))
factors_MCIA = readRDS (paste('03_Results/', set, '/01_MCIA/', o, '/factors_mcia.RDS', sep = ''))

## Extract Components 
compo_mcia = lapply(  MCIA_dry$coa , function (x) {  df = x$co ; colnames(df) = paste0("MCIA_Comp_" , 1:10) ; df } )
saveRDS (compo_mcia , file= paste('03_Results/', set, '/01_MCIA/', o, '/compo_mcia.RDS', sep = ''))

## get MCIA features 
features_w_mcia = lapply (MCIA_dry$coa , function (x) { df = x$li ; colnames(df) = paste0("MCIA_" , 1:10) ; df  })
saveRDS(features_w_mcia , file = paste('03_Results/', set, '/01_MCIA/', o, '/features_w_mcia.RDS', sep = ''))

# get MCIA compo 
compo_mcia = lapply (MCIA_dry$coa , function (x) { df = x$c1 ; colnames(df) = paste0("MCIA_Compo_" , 1:10) ; df } )

################
###   MOFA   ###
################
MOFA_dry = readRDS (paste('03_Results/', set, '/02_MOFA/', o, '/02_MOFAobject_trained_dry.RDS', sep = ''))
factors_MOFA =  readRDS (paste('03_Results/', set, '/02_MOFA/', o, '/factors_s_mofa.RDS', sep = ''))

## get MOFA  features 
features_w_mofa = get_weights (MOFA_dry)
features_w_mofa = lapply (features_w_mofa , function(x) {colnames (x) = paste0 ("MOFA_" , 1:10) ; x } )
saveRDS(features_w_mofa , file = paste('03_Results/', set, '/02_MOFA/', o, '/features_w_mofa.RDS', sep = ''))

factors_mofa = get_factors (MOFA_dry) $group1 
colnames (factors_mofa) = paste0("MOFA_" , 1:10)

## get MOFA compo 
compo_mofa = MOFA_dry@expectations$W
compo_mofa = lapply (compo_mofa , function (x) {colnames(x) = paste0("MOFA_Compo" , 1:10) ; x } )

##################
###   IntNMF   ###
##################
INTNMF_dry = readRDS (paste('03_Results/', set, '/03_IntNMF/', o, '/03_factorizations_intnmf_dry.RDS', sep = ''))

### Extract intNMF Results 

factors_intnmf = INTNMF_dry$W
colnames(factors_intnmf) = paste0( "IntNMF_" , 1:10 )

features_w_intNMF =list()

for(j in 1:length(omics [[o]])){
  features_w_intNMF[[j]] = t(INTNMF_dry$H$H1[, j]) 
  colnames(features_w_intNMF [[j]] ) = paste0( "IntNMF_" , 1:10 )
}
names (features_w_intNMF) = c("miRNA",  "mRNA" , "DNAm")

saveRDS(factors_intnmf    , file = paste('03_Results/', set, '/03_IntNMF/', o, '/factors_intnmf.RDS', sep = ''))
saveRDS(features_w_intNMF , file = paste('03_Results/', set, '/03_IntNMF/', o, '/features_w_intNMF.RDS', sep = ''))

################
###   JIVE   ###
################
JIVE_dry = readRDS (paste('03_Results/', set, '/04_JIVE/', o, '/04_factorizations_jive_dry.RDS', sep = ''))

### Extract Jive Results 
rankJV   = JIVE_dry$rankJ
rankIV.v = JIVE_dry$rankA   

J  = JIVE_dry$joint[[1]][1,]
ng = dim(JIVE_dry$joint[[1]])[1]

factors_jive = features_w_jive = list()

for(j in 2:length(omics[[o]])){
  J  = rbind(J, JIVE_dry$joint[[1]][j,]);
  ng = c(ng, dim(JIVE_dry$joint[[1]])[1])
}

svd.o = svd(J)
jV    = svd.o$v %*% diag(svd.o$d);
features_w_jive  =  svd.o$u[, 1:rankJV]
rownames(features_w_jive) = colnames( omics[[o]])
colnames(features_w_jive) = paste0("JIVE_" , 1:num_factors )

factors_jive = jV[,1:rankJV]


colnames(factors_jive) = paste0("JIVE_" , 1:10 )
rownames(factors_jive) = rownames(omics$mRNA)
names (features_w_jive) = c("miRNA",  "mRNA" , "DNAm")

saveRDS(factors_jive , file = paste('03_Results/', set, '/04_JIVE/', o, '/factors_jive.RDS', sep = ''))
saveRDS(features_w_jive , file = paste('03_Results/', set, '/04_JIVE/', o, '/features_w_jive.RDS', sep = ''))

## Extract Components JIVE

get_jiv_omic_comp = function (JIVE_dry , x) {
  svd_s = svd (JIVE_dry$joint[[x]])
  s_v = svd_s$v %*% diag(svd_s$d) 
  s = s_v [,1:10]
  colnames (s) = paste0("JIVE_Compo_" , 1:10) 
  rownames (s) = omics$covariates %>% rownames 
  return(s)
}

compo_jive = lapply (1:3 , function (x) {get_jiv_omic_comp (JIVE_dry , x )} ) 
names(compo_jive) = c("miRNA","mRNA","DNAm")
saveRDS(compo_jive , file = paste('03_Results/', set, '/04_JIVE/', o, '/compo_jive.RDS', sep = ''))

##################
###   SCIKIT   ###
##################
SCKIT_dry = readRDS (paste('03_Results/', set, '/05_SCIKIT/', o, '/05_factorizations_scikit_dry.RDS', sep = ''))

### Extract sckit Results 

factors_sckit = SCKIT_dry[[1]]
features_w_sckit = SCKIT_dry[[2]]
colnames (factors_sckit) = paste0 ( "Scikit_", 1:10)
features_w_sckit               = lapply(features_w_sckit , function (x) { colnames(x) =  paste0 ( "Scikit_", 1:10) ; x })
saveRDS(factors_sckit , file    = paste('03_Results/', set, '/05_SCIKIT/', o, '/factors_sckit.RDS', sep = ''))
saveRDS(features_w_sckit , file = paste('03_Results/', set, '/05_SCIKIT/', o, '/features_w_sckit.RDS', sep = ''))

#################
###   RGCCA   ###
#################
RGCCA_dry  = readRDS (paste('03_Results/', set, '/06_RGCCA/', o, '/06_factorizations_SGCCA_dry.RDS', sep = ''))

### Extract RGCCA Results from 4th block as matrix 

factors_RGCCA = RGCCA_dry$Y$superblock
colnames (factors_RGCCA) = paste0 ( "RGCCA_", 1:10)
saveRDS(factors_RGCCA , file = paste('03_Results/', set, '/06_RGCCA/', o, '/factors_RGCCA.RDS', sep = ''))

features_w_rgcca = list()
df_super_rgcca = as.data.frame(RGCCA_dry$a$superblock)
colnames(df_super_rgcca) = paste0("RGCCA_" , 1:10)
features_w_rgcca$miRNA =  df_super_rgcca [grep("^h",rownames(df_super_rgcca)) ,]
features_w_rgcca$mRNA =  df_super_rgcca [grep("^E",rownames(df_super_rgcca)) ,]
features_w_rgcca$DNAm =  df_super_rgcca [grep("^c",rownames(df_super_rgcca)) ,]

saveRDS(features_w_rgcca , file = paste('03_Results/', set, '/06_RGCCA/', o, '/features_w_rgcca.RDS', sep = ''))



features_w_rgcca_singleomic = RGCCA_dry$a [1:3]
features_w_rgcca_singleomic = lapply (features_w_rgcca_singleomic , function(x) {colnames (x) =  paste0("RGCCA_" , 1:10) ; x }  ) 
saveRDS(features_w_rgcca_singleomic , file = paste('03_Results/', set, '/06_RGCCA/', o, '/features_w_rgcca_singleomic.RDS', sep = ''))


### get RGCCA Compo 

compo_rgcca = RGCCA_dry$Y 
featu = lapply (compo_rgcca , function (x) { colnames(x) = paste0("RGCCA_Compo_" , 1:10) ; x  })

### get RGCCA Compo 

compo_rgcca = RGCCA_dry$Y 
featu = lapply (compo_rgcca , function (x) { colnames(x) = paste0("RGCCA_Compo_" , 1:10) ; x  })




#############################################################################################################################################






















### get the ranking of the feature based on the most correlated feature to the phenotype

## get best corr to phenotype 

cor_fact_pheno = function (compo) {
    
    p_val_cor = cor.mtest (compo)
    pval_group = p_val_cor [ 1:nrow(p_val_cor)-1, ncol(p_val_cor)] 
    
    name_fact_ = names(which.min(pval_group [ pval_group < 0.05]))
    if (length(name_fact_) != 0 ) {
        return (gsub ( "_Compo_" , "_", name_fact_ ) )
    }
    
}

Group =   as.double(as.factor(omics$covariates$Group))
all_facts = cbind(factors_intnmf,factors_jive,factors_MCIA,factors_MOFA,factors_sckit,factors_RGCCA , Group)

### Calculate an optimal ranikng based on the diffetent ones 
compos = list (factors_intnmf = factors_intnmf,
               factors_jive   = factors_jive,
               factors_MCIA   = factors_MCIA,
               factors_MOFA   = factors_MOFA,
               factors_sckit  = factors_sckit,
               factors_RGCCA  = factors_RGCCA)
compos_best =   lapply (compos , function (Y) {  cor_fact_pheno ( cbind (Y , Group) )  } )

### Result correlated according to the case control status

compos_best

gen_omic_n = function (features_w , best_fact , o , omics , frac = 0.01) {
    
    if ( !is.null (best_fact) ) { # get row names of the omic
        feat_n =  features_w[[ o ]]     %>% 
                    as.data.frame                          %>% 
                    arrange ( - abs( . [ , best_fact ]))  %>% 
                    top_frac (  frac )                     %>% 
                    row.names

        return( omics [[ o]] [ , feat_n ]  )
      }
}




fracts = list (miRNA = 0.01 , mRNA = 0.1 , DNAm =  0.1)
om = list (miRNA =  "miRNA" , mRNA =  "mRNA" , DNAm = "DNAm")

data_intNMF = mapply ( function (x , y) {  gen_omic_n (features_w_intNMF , compos_best$factors_intnmf , x , omics = omics , frac = y )}  , om , fracts  )
data_rgcca = mapply (  function (x , y) {  gen_omic_n (features_w_rgcca_singleomic , compos_best$factors_RGCCA , x , omics = omics , frac = y )}  , om , fracts  )
data_jive = mapply (   function (x , y) {  gen_omic_n (features_w_jive , compos_best$factors_jive , x , omics = omics , frac = y )}  , om , fracts  )
data_mcia = mapply (   function (x , y) {  gen_omic_n (features_w_mcia , compos_best$factors_MCIA , x , omics = omics , frac = y )}  , om , fracts  )
data_scikit = mapply ( function (x , y) {  gen_omic_n (features_w_sckit , compos_best$factors_sckit , x , omics = omics , frac = y )}  , om , fracts  )
data_mofa = mapply ( function (x , y) {  gen_omic_n (features_w_mofa , compos_best$factors_MOFA , x , omics = omics , frac = y )}  , om , fracts  )



# SNF
library(SNFtool)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(matrixStats)
library(mclust)
source(file.path (libClustpackage,'NEMO.R') )



nemo_f = function (data, tl2 = tl , names = c("miRNA", "mRNA", "DNAm", "AllOmics")) 
{
    setClass(Class = "NemoRes", representation(graph = "matrix", 
        clustering = "integer"))
    nemo.clust = function(omics.list, num.clusters = NULL, num.neighbors = NA) {
        if (is.null(num.clusters)) {
            num.clusters = NA
        }
        graph = nemo.affinity.graph(omics.list, k = num.neighbors)
        if (is.na(num.clusters)) {
            num.clusters = nemo.num.clusters(graph)
        }
        clustering = spectralClustering(graph, num.clusters)
        names(clustering) = colnames(graph)
        return(new("NemoRes", graph = graph, clustering = clustering))
    }
    nemo.clust_momix = lapply(data, function(x) nemo.clust(list(x), 
        num.clusters = 2, num.neighbors = NA))
    nemo_all = nemo.clust(data, num.clusters = 2, num.neighbors = NA)
    nemo_group = do.call(cbind, lapply(nemo.clust_momix, function(x) {
        x@clustering
    })) %>% as.data.frame
    nemo_group$ALL = nemo_all@clustering
    nemo_group$TrueLabel = tl2
    ARIs_nemo = lapply(nemo_group[, 1:(ncol(nemo_group) - 1)], 
        function(x) adjustedRandIndex(x, nemo_group[, ncol(nemo_group)]))
    names(ARIs_nemo) = names
     return(list(ARIs = ARIs_nemo, Clusters = nemo_group))
}
                       

SNF_f = function (Data,
                K = 30 ,
                alpha = 0.6 ,
                T = 50  ,
                C = 2  ,
                tl2 = tl , 
                color = c("#ffe935","#e289bf"),
                names = c("miRNA", "mRNA", "DNAm", "AllOmics")) {
    
    
    Dist.MOMIX_metaG = lapply ( lapply(Data  , t) , function (x) {(dist2((as.matrix(x)),as.matrix(x)))^(1/2)})
    W.MOMIX_metaG    = lapply (Dist.MOMIX_metaG , function (x) {affinityMatrix(x, K, alpha)})
    W_all_MOMIX_metaG           = SNF(W.MOMIX_metaG, K, T)
    group_mom             = lapply(  W.MOMIX_metaG   , function(x) {spectralClustering(x,C)} )
    group_mom             = append (group_mom , list(factALL = spectralClustering ( W_all_MOMIX_metaG  , C) ) )
    
    # get the True clusturing 
    levels(tl2)  = c(1,2)
    ##### Get Clusters 
    M_label_ = t(do.call(rbind, group_mom)) %>% as.data.frame
    M_label_$TrueLabel = tl2
    rownames(M_label_) = rownames(W.MOMIX_metaG[[1]])
    
    ARIs = lapply (  M_label_ [ , 1: (ncol(M_label_)-1 )]  , function(x) adjustedRandIndex(x , M_label_[ ,ncol(M_label_)]))
    names (ARIs) = names 
    return(  list(ARIs = ARIs ,
             Clusters = M_label_ ))
}



data_sing_omic = list ( mofa = lapply ( data_mofa ,t )  ,
                       rgcca = lapply ( data_rgcca,t )  ,
                       jive = lapply ( data_jive,t )  , 
                       mcia = lapply ( data_mcia ,t ) ,
                     #  scikit = lapply ( data_scikit,t )  ,
                       intNMF = lapply ( data_intNMF,t ) )



tl = omics$covariates$Group %>% as.factor %>% as.integer



### apply clusturing on data 

SNF_all  = lapply (  data_sing_omic , function(x) { SNF_f ( plyr::compact ( x)   , names = c( names (plyr::compact (x))  , "All"  ))  } )
NEMO_all = lapply (  data_sing_omic , function(x) { nemo_f ( plyr::compact (x)   , names = c( names (plyr::compact (x))  , "All"  ))  } )

SNF_miRNA_mRNA  = lapply (  data_sing_omic , function(x) { SNF_f (  plyr::compact (x [1:2])   , names = c( names (plyr::compact  (x[1:2]))  , "All"  ))  } )
NEMO_miRNA_mRNA = lapply (  data_sing_omic , function(x) { nemo_f (  plyr::compact (x[1:2] )   , names = c( names (plyr::compact (x[1:2]))  , "All"  ))  } )

SNF_miRNA_DNAm  = lapply (  data_sing_omic , function(x) { SNF_f (  plyr::compact (x [ c( 1,3 )])   , names = c( names (plyr::compact  (x[ c( 1,3 )]))  , "All"  ))  } )
NEMO_miRNA_DNAm = lapply (  data_sing_omic , function(x) { nemo_f (  plyr::compact (x[ c( 1,3 )] )   , names = c( names (plyr::compact (x[ c( 1,3 )]))  , "All"  ))  } )

SNF_mRNA_DNAm  = lapply (  data_sing_omic , function(x) { SNF_f (  plyr::compact (x [ c( 2,3 )])   , names = c( names (plyr::compact  (x[ c( 2,3 )]))  , "All"  ))  } )
NEMO_mRNA_DNAm = lapply (  data_sing_omic , function(x) { nemo_f ( plyr::compact (x[ c( 2,3 )] )   , names = c( names (plyr::compact (x[ c( 2,3 )]))  , "All"  ))  } )



allARIs_Supervised= plyr::rbind.fill ( lapply (  SNF_all , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  NEMO_all , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  SNF_miRNA_mRNA , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  NEMO_miRNA_mRNA , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  SNF_miRNA_DNAm , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  NEMO_miRNA_DNAm , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  SNF_mRNA_DNAm , function (x) x$ARIs ) %>% as.data.frame ,
        lapply (  NEMO_mRNA_DNAm , function (x) x$ARIs ) %>% as.data.frame  )
                
allARIs_Supervised $ method = c("SNF_all" ,
                                "NEMO_all" ,
                                "SNF_miRNA_mRNA" ,
                                "NEMO_miRNA_mRNA" ,
                                "SNF_miRNA_DNAm" ,
                                "NEMO_miRNA_DNAm" ,
                                "SNF_mRNA_DNAm" ,
                                "NEMO_mRNA_DNAm")



allARIs_Supervised = allARIs_Supervised %>% gather ( key,value , - method) 
allARIs_Supervised = cbind (allARIs_Supervised , stringr::str_split_fixed ( allARIs_Supervised$key , "[.]" , 2) )



colnames (allARIs_Supervised) =  c( "Clusturing","key","value","method","omics")



head(allARIs_Supervised)



options(repr.plot.width=20, repr.plot.height=25)

allARIs_Supervised %>% ggplot (. , aes(x=omics, y=value )) + 
    geom_boxplot( alpha = 0.2) + geom_jitter( aes(shape =  Clusturing , color = key) ,
              position=position_jitter(0.4) ,
              size = 5)+
   scale_shape_manual(values= c(16,17,18,21,22,23,24,25) , name = "Feature Selection Methods") +
    labs(x="Omics", y="ARIs", colour="Clusturing Methods ") +
coord_flip() + facet_wrap( ~ method , ncol = 3 , nrow = 7) + 
    theme(strip.text = element_text(size=18),
         axis.text=element_text(size=14),
         legend.key.size = unit(2, 'cm'),
         legend.text = element_text(size=14))


