setwd('D:/CC/Scripts/R/gitHub_MDD/multi_omics')
# A la fin: session info

library(omicade4)
library(r.jive)
library(RGCCA)
library(MOFA2)
library(IntNMF)
library(dplyr)
library(ggplot2)
library(corrplot)

################
###   DATA   ###
################

omics = readRDS ("data/MDD_Male_DNAm_No_preserv.RDS")
summary(omics)

 # Number of Factors
num.factors = 10
omics$covariates$Group = as.numeric (as.factor (omics$covariates$Group))

################
###   MCIA   ###
################

factorizations_mcia = mcia( lapply(omics [1:3] ,t)  , cia.nf = num.factors)

### Get the factors Matrix
factors_mcia = as.matrix(factorizations_mcia$mcoa$SynVar)
colnames (factors_mcia) = paste0( "MCIA_", 1:10)

compo_mcia = lapply(  factorizations_mcia$coa , function (x) {  df = x$co ; colnames(df) = paste0("MCIA_Comp_" , 1:10) ; df } )

features_w_mcia = lapply (factorizations_mcia$coa , function (x) { df = x$li ; colnames(df) = paste0("MCIA_" , 1:10) ; df  })

factors_mcia_cov = cbind(factors_mcia, omics$covariates[ match(rownames(factors_mcia) , rownames( omics$covariates) ) ,] )

rm(factorizations_mcia)

################
###   JIVE   ###
################
factorizations_jive = jive( lapply (omics[1:3], t) , 
                            rankJ   = num.factors, 
                            rankA   = rep ( num.factors, length( omics ) ), 
                            method  = "given", 
                            conv    = "default", 
                            maxiter = 1000, 
                            showProgress = T)

rankJV = factorizations_jive$rankJ
rankIV.v = factorizations_jive$rankA

J  = numeric(0)
ng = 0

factors_jive = list()
features_w_jive = list()

for(j in 1:length(omics [1:3])){
  J = rbind(J, factorizations_jive$joint[[j]])
  ng = c(ng,dim(factorizations_jive$joint[[j]])[1])
}
svd.o = svd(J)
jV = svd.o$v %*% diag(svd.o$d)

for(j in 1:length(omics[1:3])){
  features_w_jive[[j]]  =  svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]
  rownames(features_w_jive[[j]]) = rownames( t(omics[[j]]))
  colnames(features_w_jive[[j]]) = paste0("JIVE_" , 1:num.factors )
}
factors_jive = jV[, 1:rankJV]

colnames(factors_jive) = paste0("JIVE_" , 1:10 )
rownames(factors_jive) = rownames(omics$mRNA)
names (features_w_jive) = c("miRNA",  "mRNA" , "DNAm")

get_jiv_omic_comp = function (factorizations_jive , x) {
  svd_s = svd (factorizations_jive$joint[[x]])
  s_v = svd_s$v %*% diag(svd_s$d)
  s = s_v [,1:10]
  colnames (s) = paste0("JIVE_Compo_" , 1:10)
  rownames (s) = omics$covariates %>% rownames
  return(s)
}
compo_jive = lapply (1:3 , function (x) {get_jiv_omic_comp (factorizations_jive , x )} ) 
names(compo_jive) = c("miRNA","mRNA","DNAm")

rm(factorizations_jive)
rm(J)
rm(jV)
rm(svd.o)
rm(j)
rm(ng)
rm(rankIV.v)
rm(rankJV)

#################
###   RGCCA   ###
#################
result_sgcca = rgcca(blocks  = lapply (omics[1:3], t)  , 
                     ncomp = rep(num.factors, length(lapply (omics[1:3], t))), 
                     scheme = "centroid", 
                     scale = TRUE, 
                     init = "svd",
                     bias = TRUE, 
                     tol = 1e-08, 
                     verbose = T,
                     superblock = T)
factors_RGCCA = result_sgcca$Y$superblock[1:58,]
colnames (factors_RGCCA) = paste0 ( "RGCCA_", 1:10)

compo_rgcca = result_sgcca$Y[1:3]
features_w_rgcca = lapply (compo_rgcca , function (x) { colnames(x) = paste0("RGCCA_Compo_" , 1:10) ; x  })

rm(result_sgcca)

##################
###   IntNMF   ###
##################
omics_pos<-list()
for(j in 1:3){
  if(min(omics[[j]])<0){
    omics_pos[[j]]<-omics[[j]]+abs(min(omics[[j]]))
  }else{
    omics_pos[[j]]<-omics[[j]]
  }
  omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
}

factorizations_intnmf = nmf.mnnals(dat = lapply (omics_pos, as.matrix ) ,  k  = num.factors)

factors_intnmf = factorizations_intnmf$W
colnames(factors_intnmf) = paste0( "IntNMF_" , 1:10 )

features_w_intNMF =list()

for(j in 1:length(omics [1:3])){
  features_w_intNMF[[j]] = t(factorizations_intnmf$H[[j]]) 
  # rownames(features_w_intNMF [[j]] ) = rownames(omics[[j]])
  colnames(features_w_intNMF [[j]] ) = paste0( "IntNMF_" , 1:10 )
}
names (features_w_intNMF) = c("miRNA",  "mRNA" , "DNAm")

como_intNMF_miRNA = readRDS ('comp_intNMFmiRNA.RDS')
como_intNMF_mRNA = readRDS ('comp_intNMFmRNA.RDS')
como_intNMF_DNAm = readRDS ('comp_intNMFDNAm.RDS')

como_intNMF = list ( miRNA = como_intNMF_miRNA$W , mRNA = como_intNMF_mRNA$W , DNAm = como_intNMF_DNAm$W )
compo_intNMF = lapply (como_intNMF , function (x) { colnames(x) = paste0("IntNMF_Compo_" , 1:10) ; x  })

rm(factorizations_intnmf)
rm(omics_pos)

################
###   MOFA   ###
################

# Creating MOFA object from a list of matrice
MOFAobject             = create_mofa( lapply (omics , t ) )
DataOptions            = get_default_data_options(MOFAobject)
model_opts             = get_default_model_options(MOFAobject)
model_opts$num_factors = num_factors 
TrainOptions           = get_default_training_options(MOFAobject)


# Prepare the MOFA Object 
MOFAobject =  prepare_mofa( object = MOFAobject,
                            data_options = DataOptions,
                            model_options = model_opts,
                            training_options = TrainOptions )

# Configure python for MOFA2
use_basilisk = FALSE
reticulate::use_python("/opt/miniconda3/bin/python3.7", required=TRUE)
outfile = file.path(getwd(),"model2.hdf5")

MOFAobject_trained = run_mofa(MOFAobject, outfile)

###Results availiable as MOFAobject_trained.rds
# MOFAobject_trained = readRDS ("02_MOFAobject_trained_dry.rds")

features_w_mofa = get_weights(MOFAobject_trained)
features_w_mofa = lapply (features_w_mofa , function(x) {colnames (x) = paste0 ("MOFA_" , 1:10) ; x } )

factors_mofa = get_factors (MOFAobject_trained) $group1 
colnames (factors_mofa) = paste0("MOFA_" , 1:10)

compo_mofa = MOFAobject_trained@expectations$W
compo_mofa = lapply (compo_mofa , function (x) {colnames(x) = paste0("MOFA_Compo" , 1:10) ; x } )

rm(MOFAobject_trained)

#########################
###   Best features   ###
#########################

Group = as.double(as.factor(omics$covariates$Group))
all_facts = cbind(factors_mcia, factors_jive, factors_RGCCA, factors_intnmf, factors_mofa, Group)

cor.mtest = function(mat, ...) {
  mat = as.matrix(mat)
  n = ncol(mat)
  p.mat = matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp = cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] = p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) = rownames(p.mat) <- colnames(mat)
  p.mat
}

M = cor(all_facts)
p.mat = cor.mtest(all_facts)


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

### Calculate an optimal ranikng based on the diffetent ones
#### Get the best factors by omics and method 
compos = list(compo_intNMF = compo_intNMF,
              compo_rgcca = compo_rgcca,
              compo_jive = compo_jive,
              compo_mcia = compo_mcia)

compos_best = lapply (compos , function (Y) { lapply ( Y , function (x) { cor_fact_pheno ( cbind (x , Group) )  } ) } )

### Result correlated according to the case control status
compos_best
### Result correlated according to the HDRS score 
compos_best_HDRS = lapply (compos , function (Y) { lapply ( Y , function (x) { cor_fact_pheno ( cbind (x , omicsHDRS_total ) )  } ) } )
compos_best_HDRS

### get the omics featurs by method
gen_omic_n = function (features_w , best_fact , o , omics , frac = 0.01) {
  
  if ( !is.null (best_fact) ) { # get row names of the omic
    feat_n =  features_w[[ o ]] %>% 
      as.data.frame %>% 
      arrange ( - abs( . [ , best_fact ])) %>% 
      top_frac (  frac ) %>% 
      row.names
      return(omics[[ o]][, feat_n] )
  }
}
fracts = list (miRNA = 0.01 , mRNA = 0.1 , DNAm =  0.1)
om = list (miRNA =  "miRNA" , mRNA =  "mRNA" , DNAm = "DNAm")

data_intNMF = mapply ( function (x , y) {  gen_omic_n (features_w_intNMF , compos_best$compo_intNMF [[x]] , x , omics = omics , frac = y )}  , om , fracts  )
data_rgcca = mapply ( function (x , y) {  gen_omic_n (features_w_rgcca_singleomic , compos_best$compo_rgcca [[x]] , x , omics = omics , frac = y )}  , om , fracts  )
data_jive = mapply ( function (x , y) {  gen_omic_n (features_w_jive , compos_best$compo_jive [[x]] , x , omics = omics , frac = y )}  , om , fracts  )
data_mcia = mapply ( function (x , y) {  gen_omic_n (features_w_mcia , compos_best$compo_mcia [[x]] , x , omics = omics , frac = y )}  , om , fracts  )

data_sing_omic = list (intNMF = data_intNMF , rgcca = data_rgcca , jive = data_jive , mcia = data_mcia )

SNF_all  = lapply (  data_sing_omic , function(x) { SNF_f ( lapply ( plyr::compact (x) , t)  , names = c( names (plyr::compact (x))  , "All"  ))  } )
NEMO_all = lapply (  data_sing_omic , function(x) { nemo_f ( lapply ( plyr::compact (x) , t)  , names = c( names (plyr::compact (x))  , "All"  ))  } )
