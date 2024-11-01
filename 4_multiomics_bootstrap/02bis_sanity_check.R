rm(list = ls())

library(dplyr)
library(RGCCA)
library(r.jive)
library(pbmcapply)

args            = commandArgs(trailingOnly=TRUE)
experiment_name = args[1]
omics           = c("mRNA", "miRNA", "DNAm")
chunks          = c(1:24, "reference")

cohorts    = c("_female_", "_male_", "_pooled_")
cohort     = names(which(sapply(cohorts, function(x) grepl(pattern = x, x = experiment_name))))
cohort     = gsub(pattern = "_", replacement = "", x = cohort)

path_corrected_files = paste("results/02_Preprocessing", cohort, sep = '/')
for (chunk in chunks){
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
  if (chunk == "reference"){
    omics_corrected = lapply(omics, function(omic) list(readRDS(corrected_files[[omic]])))
  } else {
    omics_corrected = lapply(omics, function(omic) readRDS(corrected_files[[omic]]))
  }
  names(omics_corrected) = omics
  if (length(unique(sapply(omics_corrected, length))) != 1){
    print(sapply(omics_corrected, length))
    stop("The number of file per omic is different")
  }
  for (idx in 1:length(omics_corrected[[1]])){
    if (length(unlist(apply(sapply(omics_corrected, function(x) rownames(x[[idx]])), 1, unique))) != dim(omics_corrected[[1]][[1]])[1]){
      print(apply(sapply(omics_corrected, function(x) rownames(x[[idx]])), 1, unique))
      stop("Subjects do not match accross omic blocks")
    }
  }
  
  check_fail   = sapply(omics_corrected, function(x) sapply(x, length) == 0)
  if (chunk == "reference"){
    sample_to_rm = which(!all(!check_fail))
  } else {
    sample_to_rm = which(!apply(!check_fail, 1, all))
  }
  
  if (length(sample_to_rm) !=0){
    print(sample_to_rm)
    for (omic in omics){
      for (idx in sample_to_rm){
        omics_corrected[[omic]][[idx]] = NULL
      }
      saveRDS(object = omics_corrected[[omic]], file = corrected_files[[omic]])
    }
  }
  print(chunk)
}
