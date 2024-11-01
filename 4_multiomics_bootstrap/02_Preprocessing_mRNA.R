rm(list = ls())

library(dplyr)
library(stringr)
library(DESeq2)
library(ChAMP)
library(sva)
library(pbmcapply)
library(neuroCombat)

args                   = commandArgs(trailingOnly=TRUE)
bootstrap_samples_file = args[1]
mc.cores               = as.numeric(args[2])
cohort                 = basename(dirname(bootstrap_samples_file))
experiment_name        = basename(bootstrap_samples_file)
experiment_name        = gsub(pattern     = "^bootstrap_samples_|.RDS$", 
                              replacement = "", 
                              x           = experiment_name)

### input: covariables
covariates = readRDS(file = paste0("data/", cohort, "/covariates.RDS"))

### input: cv folds
bootstrap_samples = readRDS(file = bootstrap_samples_file)

################
###   mRNA   ###
################

data = readRDS(file = paste0("data/", cohort, "/data_mRNA.RDS"))

correct_mRNA = function (data , sample , covariates, All = TRUE) {
  
  ### filter low expressed genes
  train_mRNA_f = data [rowMeans (data) > 10 ,]
  train_mRNA_f = train_mRNA_f [ , sample]
  
  ### vst norm 
  dds_train_adj = DESeqDataSetFromMatrix( (train_mRNA_f) , colData=covariates[sample,] , design=~1)
  vst_train_adj = vst(dds_train_adj)
  
  ### Correct data 
  if (All == TRUE){
    covmRNA  = c("BMI",  "Rin", "AGE", "SEX" ,"Polynuclear_neutrophile", "Lymphocyte")
    
    corrected_mRNA_train_model = (lm(t(assay(vst_train_adj))~ Rin + BMI + AGE + SEX + Polynuclear_neutrophile + Lymphocyte, 
                                     data=covariates[sample, covmRNA]) )
    
    corrected_mRNA_train = corrected_mRNA_train_model$residuals
  }else{
    covmRNA  = c("BMI",  "Rin", "AGE" ,"Polynuclear_neutrophile", "Lymphocyte")
    
    corrected_mRNA_train_model = (lm(t(assay(vst_train_adj))~ Rin + BMI + AGE  + Polynuclear_neutrophile + Lymphocyte, 
                                     data=covariates[sample, covmRNA]) )
    
    corrected_mRNA_train = corrected_mRNA_train_model$residuals
  }
  return (corrected_mRNA_train)
}

#For correcting all folds:
if (cohort == "pooled"){
  bootstrap_mRNA_corrected = pbmclapply(bootstrap_samples, function(x) correct_mRNA(data, x, covariates), mc.cores = mc.cores)
} else {
  bootstrap_mRNA_corrected = pbmclapply(bootstrap_samples, function(x) correct_mRNA(data, x, covariates, FALSE), mc.cores = mc.cores)
}

bootstrap_sample_with_nullVar = sapply(bootstrap_mRNA_corrected, function(x) all(apply(x, 2, sd) != 0))
if (length(which(!bootstrap_sample_with_nullVar)) != 0){
  for (idx in which(!bootstrap_sample_with_nullVar)){
    bootstrap_mRNA_corrected[[idx]] = bootstrap_mRNA_corrected[[idx]]*NULL
  }
}

all_colnames = sapply(bootstrap_mRNA_corrected, colnames)
if (all(class(all_colnames) == "list")){
  print(all_colnames)
  stop("There is a different number of variables per bootstrap sample")
}

unique_colnames = apply(all_colnames, 1, unique)
if(class(unique_colnames) == "list"){
  print(unique_colnames)
  stop("Variables are not the same or in the same ordre accross bootstrap samples")
}



chunk_size = 48
B          = length(bootstrap_samples)
chunk_idx  = seq(from = 1, to = B, by = chunk_size)
for (nb_chunk in 1:length(chunk_idx)){
  chunk_start = chunk_idx[nb_chunk]
  if (nb_chunk == length(chunk_idx)){
    chunk_end = B
  } else {
    chunk_end = chunk_idx[nb_chunk + 1] - 1
  }
  saveRDS(bootstrap_mRNA_corrected[chunk_start:chunk_end], file = paste0("results/02_Preprocessing/", 
                                                                         cohort, 
                                                                         "/mRNA/bootstrap_mRNA_corrected_", 
                                                                         experiment_name,
                                                                         "_nbChunk_",
                                                                         nb_chunk,
                                                                         ".RDS"))
}

if (cohort == "pooled"){
  reference_mRNA_corrected = correct_mRNA(data, colnames(data), covariates)
} else {
  reference_mRNA_corrected = correct_mRNA(data, colnames(data), covariates, FALSE)
}

if (is.null(reference_mRNA_corrected)){
  stop("Reference sample generated NULL output")
}

if (length(reference_mRNA_corrected) == 0){
  stop("Issue with the reference sample")
}

if (sum(colnames(reference_mRNA_corrected) == colnames(bootstrap_mRNA_corrected[[1]]))/dim(reference_mRNA_corrected)[2] !=1){
  stop("There is a different number of variables between the reference smaple and all the other bootstrap samples")
}


saveRDS(reference_mRNA_corrected , file = paste0("results/02_Preprocessing/", 
                                                 cohort, 
                                                 "/mRNA/bootstrap_mRNA_corrected_", 
                                                 experiment_name,
                                                 "_nbChunk_",
                                                 "reference",
                                                 ".RDS"))

