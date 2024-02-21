library(dplyr)
library(RGCCA)
library(omicade4)
library(r.jive)
library(MOFA2)
library(IntNMF)

cv_DNAm_corr  = readRDS("results/2_PreProcessing/cv_DNAm_corr.RDS")
cv_miRNA_corr = readRDS("results/2_PreProcessing/cv_miRNA_corr.RDS")
cv_mRNA_corr  = readRDS("results/2_PreProcessing/cv_mRNA_corr.RDS")
cov_pooled    = readRDS(file = "data/cov_pooled.RDS")

gen_omic_n = function (features_w , best_fact , o , omics , frac = 0.1) {
  if ( !is.null (best_fact) ) { # get row names of the omic
    feat_n = features_w[[ o ]][ , best_fact ]
    feat_n = feat_n %>% as.data.frame
    feat_n = feat_n %>% abs %>% top_frac (  frac )
    feat_n = feat_n %>% row.names
    return( omics [[ o]] [ , feat_n ]  )
  }
}

data_train = data_test = list()
for (i in 1:length(cv_DNAm_corr)){
  omics_train = list(miRNA      = cv_miRNA_corr[[i]]$corrected_miRNA_train,
                     mRNA       = cv_mRNA_corr[[i]]$corrected_mRNA_train,
                     DNAm       = cv_DNAm_corr[[i]]$corrected_DNAm_train,
                     covariates = cov_pooled[rownames(cv_miRNA_corr[[i]]$corrected_miRNA_train), paste0("Hamilton_scale_item_", 1:17)])
  
  omics_test = list(miRNA      = cv_miRNA_corr[[i]]$corrected_miRNA_test,
                    mRNA       = cv_mRNA_corr[[i]]$corrected_mRNA_test,
                    DNAm       = cv_DNAm_corr[[i]]$corrected_DNAm_test,
                    covariates = cov_pooled[rownames(cv_miRNA_corr[[i]]$corrected_miRNA_test), paste0("Hamilton_scale_item_", 1:17)])
  
  # Number of Factors
  num.factors = 4
  
  ######################
  ###   RGCCA Start  ###
  ######################
  
  result_rgcca = rgcca(blocks  = omics_train,
                       ncomp = num.factors,
                       scheme = "centroid",
                       scale = TRUE,
                       init = "random", # for the sake of computational time this parameter was set in this example to "random", though, results were computed with "init = "svd", which takes more time.
                       bias = TRUE,
                       tol = 1e-08,
                       verbose = T,
                       superblock = T)
  
  features_w       = result_rgcca$a[1:3]
  features_w$miRNA = features_w$miRNA[colnames(omics_train$miRNA),]
  features_w$mRNA  = features_w$mRNA[colnames(omics_train$mRNA),]
  features_w$DNAm  = features_w$DNAm[colnames(omics_train$DNAm),]
  
  features_w = lapply(features_w , function (x) { colnames(x) = paste0("RGCCA_" , 1:num.factors) ; x  })
  
  factors           = result_rgcca$Y$superblock
  colnames(factors) = paste0 ( "RGCCA_", 1:num.factors)
  
  ####################
  ###   RGCCA End  ###
  ####################
  
  #####################################################
  ######  UNCOMMENT THE METHOD YOU WANT TO TEST   #####
  ######  (AND COMMENT ALL THE PART ON RGCCA)     ##### 
  #####################################################
  
  
  # ####################
  # ###   MCIA Start ###
  # ####################
  # 
  # MCIA_dry = mcia(lapply(omics_train[1:3], t)  , cia.nf = num.factors)
  # 
  # # Extract factors
  # factors = as.matrix(MCIA_dry$mcoa$SynVar)
  # colnames (factors) = paste0 ( "MCIA_", 1:num.factors)
  # 
  # # Extract features
  # features_w = lapply (MCIA_dry$coa , function (x) { df = x$li ; colnames(df) = paste0("MCIA_" , 1:num.factors) ; df  })
  # 
  # ##################
  # ###   MCIA End ###
  # ##################
  # 
  # ####################
  # ###   JIVE Start ###
  # ####################
  # 
  # factorizations_jive = jive(lapply (omics_train[1:3], t) , 
  #                                    rankJ        = num.factors, 
  #                                    rankA        = rep ( num.factors, length( omics_train ) ), 
  #                                    method       = "given", 
  #                                    conv         = "default", 
  #                                    maxiter      = 1000, 
  #                                    showProgress = T)
  # rankJV   = factorizations_jive$rankJ
  # rankIV.v = factorizations_jive$rankA
  # 
  # J  = numeric(0)
  # ng = 0
  # 
  # factors = list()
  # features_w = list()
  # 
  # for(j in 1:length(omics_train [1:3])){
  #   J = rbind(J, factorizations_jive$joint[[j]])
  #   ng = c(ng,dim(factorizations_jive$joint[[j]])[1])
  # }
  # svd.o = svd(J)
  # jV = svd.o$v %*% diag(svd.o$d)
  # 
  # for(j in 1:length(omics_train[1:3])){
  #   features_w[[j]]  =  svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]
  #   rownames(features_w[[j]]) = rownames( t(omics_train[[j]]))
  #   colnames(features_w[[j]]) = paste0("JIVE_" , 1:num.factors )
  # }
  # factors = jV[, 1:rankJV]
  # 
  # colnames(factors) = paste0("JIVE_" , 1:num.factors )
  # rownames(factors) = rownames(omics_train$mRNA)
  # names (features_w) = c("miRNA",  "mRNA" , "DNAm")
  # 
  # get_jiv_omic_comp = function (factorizations_jive , x) {
  #   svd_s = svd (factorizations_jive$joint[[x]])
  #   s_v = svd_s$v %*% diag(svd_s$d)
  #   s = s_v [,1:num.factors]
  #   colnames (s) = paste0("JIVE_Compo_" , 1:num.factors)
  #   rownames (s) = omics_train$covariates %>% rownames
  #   return(s)
  # }
  # compo_jive = lapply (1:3 , function (x) {get_jiv_omic_comp (factorizations_jive , x )} ) 
  # names(compo_jive) = c("miRNA","mRNA","DNAm")
  # 
  # rm(J)
  # rm(jV)
  # rm(svd.o)
  # rm(j)
  # rm(ng)
  # rm(rankIV.v)
  # rm(rankJV)
  # 
  # ##################
  # ###   JIVE End ###
  # ##################
  # 
  # ######################
  # ###   IntNMF Start ###
  # ######################
  # 
  # omics_pos<-list()
  # for(j in 1:3){
  #   if(min(omics_train[[j]])<0){
  #     omics_pos[[j]]<-omics_train[[j]]+abs(min(omics_train[[j]]))
  #   }else{
  #     omics_pos[[j]]<-omics_train[[j]]
  #   }
  #   omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
  # }
  # 
  # factorizations_intnmf = nmf.mnnals(dat = lapply (omics_pos, as.matrix ) ,  k  = num.factors)
  # 
  # factors = factorizations_intnmf$W
  # colnames(factors) = paste0( "IntNMF_" , 1:num.factors )
  # 
  # features_w =list()
  # 
  # for(j in 1:length(omics_train [1:3])){
  #   features_w[[j]] = t(factorizations_intnmf$H[[j]]) 
  #   colnames(features_w [[j]] ) = paste0( "IntNMF_" , 1:num.factors )
  # }
  # names (features_w) = c("miRNA",  "mRNA" , "DNAm")
  # 
  # rm(omics_pos)
  # 
  # ####################
  # ###   IntNMF End ###
  # ####################
  # 
  # #####################
  # ###   MOFA Start  ###
  # #####################
  # 
  # # Creating MOFA object from a list of matrice
  # MOFAobject             = create_mofa( lapply (omics_train[1:3] , t ) )
  # DataOptions            = get_default_data_options(MOFAobject)
  # model_opts             = get_default_model_options(MOFAobject)
  # model_opts$num_factors = num.factors 
  # TrainOptions           = get_default_training_options(MOFAobject)
  # 
  # # Prepare the MOFA Object 
  # MOFAobject =  prepare_mofa(object           = MOFAobject,
  #                            data_options     = DataOptions,
  #                            model_options    = model_opts,
  #                            training_options = TrainOptions)
  # 
  # # Configure python for MOFA2
  # use_basilisk = FALSE
  # reticulate::use_python("/opt/miniconda3/bin/python3.7", required=TRUE)
  # outfile = file.path(getwd(),"model2.hdf5")
  # 
  # MOFAobject_trained = run_mofa(MOFAobject, outfile)
  # 
  # features_w = get_weights(MOFAobject_trained)
  # features_w = lapply (features_w , function(x) {colnames (x) = paste0 ("MOFA_" , 1:num.factors) ; x } )
  # 
  # factors = get_factors (MOFAobject_trained) $group1 
  # colnames (factors) = paste0("MOFA_" , 1:num.factors)
  # 
  # compo_mofa = MOFAobject_trained@expectations$W
  # compo_mofa = lapply (compo_mofa , function (x) {colnames(x) = paste0("MOFA_Compo" , 1:num.factors) ; x } )
  # 
  # 
  # ###################
  # ###   MOFA End  ###
  # ###################
  
  
  
  #########################
  ###   Best Features   ###
  #########################
  
  ### get the ranking of the feature based on the most correlated feature to the phenotype
  
  ## get best corr to phenotype 
  Group    = as.double(as.factor(cov_pooled[rownames(cv_miRNA_corr[[i]]$corrected_miRNA_train), 'GROUP']))
  best_cor = colnames(factors)[which.max(abs(apply(factors, 2, function(x){cor(x, Group)})))]
  print(paste("BEST CORR :", best_cor))
  
  fracts = list (miRNA = 0.1 , mRNA = 0.1 , DNAm =  0.1)
  om     = list (miRNA =  "miRNA" , mRNA =  "mRNA" , DNAm = "DNAm")
  
  data_train[[i]] = mapply(function(x , y) {gen_omic_n (features_w , best_cor , x , omics = omics_train , frac = y )}  , om , fracts  )
  
  data_test[[i]]       = omics_test
  data_test[[i]]$mRNA  = data_test[[i]]$mRNA[, colnames(data_train[[i]]$mRNA)]
  data_test[[i]]$miRNA = data_test[[i]]$miRNA[, colnames(data_train[[i]]$miRNA)]
  data_test[[i]]$DNAm  = data_test[[i]]$DNAm[, colnames(data_train[[i]]$DNAm)]

}

saveRDS(data_train, 'results/3_FeaturesSelection/data_train.RDS')
saveRDS(data_test, 'results/3_FeaturesSelection/data_test.RDS')
