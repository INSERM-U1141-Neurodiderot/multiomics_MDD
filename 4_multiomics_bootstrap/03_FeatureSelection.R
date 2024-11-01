library(dplyr)
library(RGCCA)
library(r.jive)
library(pbmcapply)

args            = commandArgs(trailingOnly=TRUE)
experiment_name = args[1]
mc.cores        = as.numeric(args[2])
jDR_method      = args[3]
omics           = unlist(strsplit(x = args[4], split = "_"))
chunk           = args[5]

num.factors  = 10
ratioVarKept = 0.1

cohorts    = c("_female_", "_male_", "_pooled_")
cohort     = names(which(sapply(cohorts, function(x) grepl(pattern = x, x = experiment_name))))
cohort     = gsub(pattern = "_", replacement = "", x = cohort)
covariates = readRDS(file = paste0("data/", cohort, "/covariates.RDS"))

path_corrected_files = paste("results/02_Preprocessing", cohort, sep = '/')
pattern              = paste0(experiment_name, "_nbChunk_", chunk, ".RDS")
corrected_files      = c()
for (omic in  omics){
  tmp = list.files(path = paste(path_corrected_files, omic, sep = "/"), pattern = pattern, full.names = T)
  if(length(tmp) == 0){
    print(pattern)
    stop("No file found")
  }else if (length(tmp) >1) {
    print(pattern)
    stop("More than one file found")
  } else {
    corrected_files[[omic]] = tmp
  }
}

gen_omic_n = function (features_w , best_fact , o , frac = 0.1) {
  if ( !is.null (best_fact) ) { # get row names of the omic
    feat_n = features_w[[ o ]][ , best_fact ]
    feat_n = feat_n %>% as.data.frame
    feat_n = feat_n %>% abs %>% top_frac (frac)
    feat_n = feat_n %>% row.names
    return(feat_n)
  }
}

jdr = function(blocks, jDR_method, num.factors, omics){
  if (jDR_method == "RGCCA"){
    jdr_model = rgcca(blocks     = blocks,
                      ncomp      = num.factors,
                      scheme     = "centroid",
                      scale      = TRUE,
                      init       = "random", # for the sake of computational time this parameter was set in this example to "random", though, results were computed with "init = "svd", which takes more time.
                      bias       = TRUE,
                      tol        = 1e-08,
                      verbose    = T,
                      superblock = T)
    
    features_w        = jdr_model$a[omics]
    features_w        = lapply(features_w , function (x) { colnames(x) = paste0("RGCCA_" , 1:num.factors) ; x  })
    factors           = jdr_model$Y$superblock
    colnames(factors) = paste0 ( "RGCCA_", 1:num.factors)
    jdr_model$call$blocks = NULL
    jdr_model$blocks      = NULL
  } else if (jDR_method == "JIVE"){
    jdr_model = jive(data         = lapply(blocks, t) ,
                     rankJ        = num.factors,
                     rankA        = rep (num.factors, length(omics) + 1),
                     method       = "given",
                     conv         = "default",
                     maxiter      = 1000,
                     showProgress = T)
    rankJV     = jdr_model$rankJ
    J          = numeric(0)
    ng         = 0
    factors    = list()
    features_w = list()
    
    for(j in 1:(length(omics))){#+1
      J  = rbind(J, jdr_model$joint[[j]])
      ng = c(ng,dim(jdr_model$joint[[j]])[1])
    }
    
    svd.o = svd(J)
    jV    = svd.o$v %*% diag(svd.o$d)
    
    for(j in 1:length(omics)){
      features_w[[j]]           = svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]
      rownames(features_w[[j]]) = colnames(blocks[[j]])
      colnames(features_w[[j]]) = paste0("JIVE_" , 1:num.factors )
    }
    factors = jV[, 1:rankJV]
    
    colnames(factors) = paste0("JIVE_" , 1:num.factors )
    rownames(factors) = rownames(blocks[[1]])
    names (features_w) = omics
    
    rm(J)
    rm(jV)
    rm(svd.o)
    rm(j)
    rm(ng)
    rm(rankJV)
    jdr_model$data       = NULL
    jdr_model$joint      = NULL
    jdr_model$individual = NULL
  }
  duplicated_rownames  = sapply(strsplit(x = rownames(blocks[[1]]), split = "\\."), function(x) x[1])
  Group                = as.double(as.factor(covariates[duplicated_rownames, 'GROUP']))
  best_cor             = colnames(factors)[which.max(abs(apply(factors, 2, function(x){cor(x, Group)})))]
  print(paste("BEST CORR :", best_cor))
  
  selectedFeatures        = lapply(omics, function(omic) gen_omic_n(features_w = features_w, best_fact = best_cor, o = omic, frac = ratioVarKept))
  names(selectedFeatures) = omics
  return(list(jdr_model        = jdr_model, 
              features_w       = features_w, 
              factors          = factors, 
              best_cor         = best_cor, 
              selectedFeatures = selectedFeatures, 
              jDR_method       = jDR_method))
}

if (chunk == "reference"){
  omics_corrected = lapply(omics, function(omic) list(readRDS(corrected_files[[omic]])))
} else {
  omics_corrected = lapply(omics, function(omic) readRDS(corrected_files[[omic]]))
}
names(omics_corrected) = omics
blocks_per_samples = lapply(1:length(omics_corrected[[1]]), function(sample){
  cur_sample = c()
  for (omic in omics){
    tmp = omics_corrected[[omic]][[sample]]
    if (is.null(tmp)){
      return(NULL)
    }
    cur_sample[[omic]] = tmp
  }
  rownames_per_omic        = sapply(cur_sample, rownames)
  unique_rownames_per_omic = apply(rownames_per_omic, 1, unique)
  if (class(unique_rownames_per_omic) == "character"){
    print("All omics have the same individuals on the same rows")
  } else {
    print(unique_rownames_per_omic)
    stop("Individuals are not the same accross omics")
  }
  duplicated_rownames  = sapply(strsplit(x = unique_rownames_per_omic, split = "\\."), function(x) x[1])
  cur_sample[["hdrs"]] = covariates[duplicated_rownames, paste0("Hamilton_scale_item_", 1:17)]
  return(cur_sample)
})

idx_to_rm = which(sapply(blocks_per_samples, is.null))
if (length(idx_to_rm)!=0){
  print(idx_to_rm)
  stop("There shouldn't be any NULL elements")
}
rm(omics_corrected)

jdr_results = pbmclapply(blocks_per_samples, function(blocks) jdr(blocks      = blocks, 
                                                                  jDR_method  = jDR_method, 
                                                                  num.factors = num.factors, 
                                                                  omics       = omics), mc.cores = mc.cores)

saveRDS(jdr_results, file = paste0("results/03_FeatureSelection/", 
                                   cohort, 
                                   "/bootstrap_", 
                                   experiment_name, 
                                   "_jdr_",
                                   jDR_method,
                                   "_omics_",
                                   args[4],
                                   "_nbChunk_",
                                   chunk,
                                   ".RDS"))
