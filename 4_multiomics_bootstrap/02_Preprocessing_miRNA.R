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

#################
###   miRNA   ###
#################

raw.miRNA = readRDS(file = paste0("data/", cohort, "/data_miRNA.rds"))

correct_miRNA = function (data , sample , covariates, All = T) {
  
  Controls_train = covariates %>% filter ( GROUP == "Control" ) %>% rownames
  Patient_train  = covariates %>% filter ( GROUP == "MDD" ) %>% rownames
  
  # filter miRs : a mir is retained only if present with more then 1
  # reads in more then 60% of either patients or controls 
  
  x_train = rowSums(data[,Patient_train]  > 1) >= 0.6*length(Patient_train)
  y_train = rowSums(data[,Controls_train] > 1) >= 0.6*length(Controls_train)
  raw.miRNA.f.train = data [x_train | y_train , sample]
  
  ### vst norm 
  dds_train_adj = DESeqDataSetFromMatrix( (raw.miRNA.f.train) , colData=covariates[sample,] , design=~1)
  vst_train_adj = varianceStabilizingTransformation(dds_train_adj)
  
  ### Correct data
  if (All == TRUE){
    covmiRNA  = c( "BMI",  "Rin.miR", "AGE", "SEX" ,"Polynuclear_neutrophile", "Lymphocyte")
    
    corrected_miRNA_train = (lm(t(assay(vst_train_adj))~ Rin.miR + BMI + AGE + SEX + Polynuclear_neutrophile + Lymphocyte, 
                                data = covariates[sample , covmiRNA]) )$residuals
  }else{
    covmiRNA  = c( "BMI",  "Rin.miR", "AGE" ,"Polynuclear_neutrophile", "Lymphocyte")
    
    corrected_miRNA_train = (lm(t(assay(vst_train_adj))~ Rin.miR + BMI + AGE + Polynuclear_neutrophile + Lymphocyte, 
                                data = covariates[sample , covmiRNA]) )$residuals
  }
  return (corrected_miRNA_train)
}

#For correcting all folds:
if (cohort == "pooled"){
  bootstrap_miRNA_corrected = pbmclapply(bootstrap_samples, function(x) correct_miRNA(raw.miRNA, x, covariates), mc.cores = mc.cores)
} else {
  bootstrap_miRNA_corrected = pbmclapply(bootstrap_samples, function(x) correct_miRNA(raw.miRNA, x, covariates, FALSE), mc.cores = mc.cores)
}


bootstrap_sample_with_nullVar = sapply(bootstrap_miRNA_corrected, function(x) all(apply(x, 2, sd) != 0))
if (length(which(!bootstrap_sample_with_nullVar)) != 0){
  for (idx in which(!bootstrap_sample_with_nullVar)){
    bootstrap_miRNA_corrected[[idx]] = bootstrap_miRNA_corrected[[idx]]*NULL
  }
}


all_colnames = sapply(bootstrap_miRNA_corrected, colnames)
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
  saveRDS(bootstrap_miRNA_corrected[chunk_start:chunk_end] , file = paste0("results/02_Preprocessing/", 
                                                                           cohort, 
                                                                           "/miRNA/bootstrap_miRNA_corrected_", 
                                                                           experiment_name, 
                                                                           "_nbChunk_",
                                                                           nb_chunk,
                                                                           ".RDS"))
}


if (cohort == "pooled"){
  reference_miRNA_corrected = correct_miRNA(raw.miRNA, colnames(raw.miRNA), covariates)
} else {
  reference_miRNA_corrected = correct_miRNA(raw.miRNA, colnames(raw.miRNA), covariates, FALSE)
}

if (is.null(reference_miRNA_corrected)){
  stop("Reference sample generated NULL output")
}

if (length(reference_miRNA_corrected) == 0){
  stop("Issue with the reference sample")
}

if (sum(colnames(reference_miRNA_corrected) == colnames(bootstrap_miRNA_corrected[[1]]))/dim(reference_miRNA_corrected)[2] !=1){
  stop("There is a different number of variables between the reference smaple and all the other bootstrap samples")
}


saveRDS(reference_miRNA_corrected , file = paste0("results/02_Preprocessing/", 
                                                  cohort, 
                                                  "/miRNA/bootstrap_miRNA_corrected_", 
                                                  experiment_name, 
                                                  "_nbChunk_",
                                                  "reference",
                                                  ".RDS"))
