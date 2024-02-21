library(dplyr)
library(RGCCA)

cv_DNAm_corr  = readRDS("results/2_PreProcessing/cv_DNAm_corr_filter_10.RDS")
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

data_RGCCA_train = data_RGCCA_test = list()
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
  
  #################
  ###   RGCCA   ###
  #################
  result_sgcca = rgcca(blocks  = omics_train,
                       ncomp = num.factors,
                       scheme = "centroid",
                       scale = TRUE,
                       init = "random",
                       bias = TRUE,
                       tol = 1e-08,
                       verbose = T,
                       superblock = T)
  
  features_w_rgcca       = result_sgcca$a[1:3]
  features_w_rgcca$miRNA = features_w_rgcca$miRNA[colnames(omics_train$miRNA),]
  features_w_rgcca$mRNA  = features_w_rgcca$mRNA[colnames(omics_train$mRNA),]
  features_w_rgcca$DNAm  = features_w_rgcca$DNAm[colnames(omics_train$DNAm),]
  
  features_w_rgcca = lapply(features_w_rgcca , function (x) { colnames(x) = paste0("RGCCA_" , 1:num.factors) ; x  })
  
  factors_RGCCA            = result_sgcca$Y$superblock
  colnames (factors_RGCCA) = paste0 ( "RGCCA_", 1:num.factors)
  
  
  #########################
  ###   Best Features   ###
  #########################
  
  ### get the ranking of the feature based on the most correlated feature to the phenotype
  
  ## get best corr to phenotype 
  Group    = as.double(as.factor(cov_pooled[rownames(cv_miRNA_corr[[i]]$corrected_miRNA_train), 'GROUP']))
  best_cor = colnames(factors_RGCCA)[which.max(abs(apply(factors_RGCCA, 2, function(x){cor(x, Group)})))]
  print(paste("BEST CORR :", best_cor))
  
  fracts = list (miRNA = 0.1 , mRNA = 0.1 , DNAm =  0.1)
  om     = list (miRNA =  "miRNA" , mRNA =  "mRNA" , DNAm = "DNAm")
  
  data_RGCCA_train[[i]] = mapply(function(x , y) {gen_omic_n (features_w_rgcca , best_cor , x , omics = omics_train , frac = y )}  , om , fracts  )
  
  data_RGCCA_test[[i]]       = omics_test
  data_RGCCA_test[[i]]$mRNA  = data_RGCCA_test[[i]]$mRNA[, colnames(data_RGCCA_train[[i]]$mRNA)]
  data_RGCCA_test[[i]]$miRNA = data_RGCCA_test[[i]]$miRNA[, colnames(data_RGCCA_train[[i]]$miRNA)]
  data_RGCCA_test[[i]]$DNAm  = data_RGCCA_test[[i]]$DNAm[, colnames(data_RGCCA_train[[i]]$DNAm)]

}

saveRDS(data_RGCCA_train, 'results/3_FeaturesSelection/data_RGCCA_train_filter_10.RDS')
saveRDS(data_RGCCA_test, 'results/3_FeaturesSelection/data_RGCCA_test_filter_10.RDS')
