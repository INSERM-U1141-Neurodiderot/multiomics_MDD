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
mc.cores               = as.integer(args[2])
cohort                 = basename(dirname(bootstrap_samples_file))
experiment_name        = basename(bootstrap_samples_file)
experiment_name        = gsub(pattern     = "^bootstrap_samples_|.RDS$", 
                              replacement = "", 
                              x           = experiment_name)
slide                  = as.logical(args[3])
print(slide)

### input: covariables
covariates = readRDS(file = paste0("data/", cohort, "/covariates.RDS"))

### input: cv folds
bootstrap_samples = readRDS(file = bootstrap_samples_file)

################
###   DNAm   ###
################

var_filter = function (DNAm_all , freq = 0.1 )  {
  ### filter according to variance 
  rvars         = rowVars (DNAm_all)
  names(rvars)  = rownames(DNAm_all)
  top_var       = rvars %>% .[order(-.)] %>% head (n=round ( freq * length(.) ))  
  DNAm_filtered =  DNAm_all [ names(top_var) , ]
  return(DNAm_filtered)
}

neurocombat_correct_boot = function(DataTrain , CovData , Cov, boot, mean.only = F) {
  
  if (length(Cov) == 1 ) {    
    X_train = neuroCombat(dat = DataTrain, batch = as.factor(CovData[boot, Cov]), mean.only = mean.only)      
    X_train = X_train$dat.combat
  }else{
    error("Covariates of length greater than 1 is not dealt with")  
  }
  
  return( X_train )
}


correction_DNAm = function(bootstrap_samples, DNAm.npy, i, pd_mdd2, LeucocyteFraction.mdd, freq, CovDNAm, cohort, slide = F){
  Combat_DNAm = var_filter(DNAm.npy, freq)
  Combat_DNAm = Combat_DNAm[, bootstrap_samples[[i]]]
  print(i)
  for (Cov in CovDNAm) {
    if ((cohort == "female") & (Cov == "Age_bin")){
      Combat_DNAm = neurocombat_correct_boot(DataTrain  = Combat_DNAm, 
                                             CovData    = pd_mdd2,  
                                             Cov        = Cov, 
                                             boot       = bootstrap_samples[[i]], 
                                             mean.only  = T)
    } else {
      Combat_DNAm = neurocombat_correct_boot(DataTrain  = Combat_DNAm, 
                                             CovData    = pd_mdd2,  
                                             Cov        = Cov, 
                                             boot       = bootstrap_samples[[i]])
    }
  }
  
  if (slide){
    Combat_DNAm = neurocombat_correct_boot(DataTrain  = Combat_DNAm, 
                                           CovData    = pd_mdd2,  
                                           Cov        = "Slide", 
                                           boot       = bootstrap_samples[[i]],
                                           mean.only  = T)
  }
  
  DNAm_Train = as.data.frame(t(Combat_DNAm))
  
  colnames(DNAm_Train) = rownames(Combat_DNAm)
  rownames(DNAm_Train) = gsub(pattern = "^X", replacement = "", x = rownames(DNAm_Train))
  
  b.value.mdd.bin = NULL
  beta.lm = lm (as.matrix(DNAm_Train)  ~ CD4+CD8+MO+B+NK+GR , data = LeucocyteFraction.mdd[bootstrap_samples[[i]]  , 1:6 ])
  b.value.mdd.bin <- beta.lm$residuals + matrix(apply(DNAm_Train, 2, mean) ,nrow=nrow( beta.lm$residuals ), ncol=ncol( beta.lm$residuals ))
  
  b.value.mdd.bin[b.value.mdd.bin >= 1] <- 0.99999999
  b.value.mdd.bin[b.value.mdd.bin <= 0] <- 0.00000001
  
  DNAm_Train_c = b.value.mdd.bin
  
  return (DNAm_Train_c)
}

myNorm.mdd            = readRDS(file = paste0("data/", cohort, "/data_DNAm.rds"))
myNorm.mdd            = myNorm.mdd$beta
pd_mdd                = readRDS(paste0("data/", cohort, "/pd_mdd.RDS")) # pd file containes metadata of samples
LeucocyteFraction.mdd = readRDS("data/pooled/LeucocyteFraction.mdd.RDS") # leucocyte fractions estimation using Houseman method

colnames(myNorm.mdd) = rownames(pd_mdd)[match(colnames(myNorm.mdd), pd_mdd$Sample_Name)]
idx_to_rm_1          = which(is.na(colnames(myNorm.mdd)))
if (length(idx_to_rm_1) != 0){
  myNorm.mdd = myNorm.mdd[, -idx_to_rm_1]
}

tmp_rownames                    = rownames(pd_mdd)[match(rownames(LeucocyteFraction.mdd), pd_mdd$Sample_Name)]
idx_to_rm                       = which(is.na(tmp_rownames))
if (length(idx_to_rm) != 0){
  LeucocyteFraction.mdd = LeucocyteFraction.mdd[-idx_to_rm, ]
}
rownames(LeucocyteFraction.mdd) = tmp_rownames[-idx_to_rm]

if (cohort == "pooled"){
  CovDNAm  = c("Array",  "Sex", "BMI.bin", "Age_bin")
} else {
  CovDNAm = c("Array", "BMI.bin", "Age_bin")
}

bootstrap_DNAm_corrected = pbmclapply(1:length(bootstrap_samples), function(x) correction_DNAm(bootstrap_samples     = bootstrap_samples,
                                                                                               DNAm.npy              = myNorm.mdd,
                                                                                               i                     = x,
                                                                                               pd_mdd2               = pd_mdd,
                                                                                               LeucocyteFraction.mdd = LeucocyteFraction.mdd,
                                                                                               freq                  = 0.1,
                                                                                               CovDNAm               = CovDNAm,
                                                                                               cohort                = cohort,
                                                                                               slide                 = slide), mc.cores = mc.cores)

bootstrap_sample_with_nullVar = sapply(bootstrap_DNAm_corrected, function(x) all(apply(x, 2, sd) != 0))
if (length(which(!bootstrap_sample_with_nullVar)) != 0){
  for (idx in which(!bootstrap_sample_with_nullVar)){
    bootstrap_DNAm_corrected[[idx]] = bootstrap_DNAm_corrected[[idx]]*NULL
  }
}


idx_NULL     = which(sapply(bootstrap_DNAm_corrected, length) != 0)
all_colnames = sapply(bootstrap_DNAm_corrected[idx_NULL], colnames)
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
  saveRDS(bootstrap_DNAm_corrected[chunk_start:chunk_end] , file = paste0("results/02_Preprocessing/",
                                                                          cohort,
                                                                          "/DNAm/bootstrap_DNAm_corrected_",
                                                                          experiment_name,
                                                                          "_nbChunk_",
                                                                          nb_chunk,
                                                                          ".RDS"))
}

bootstrap_DNAm_corrected = readRDS(file = paste0("results/02_Preprocessing/", 
                                                  cohort,
                                                  "/DNAm/bootstrap_DNAm_corrected_",
                                                  experiment_name,
                                                  "_nbChunk_",
                                                  1,
                                                  ".RDS"))

### input: covariables
covariates = readRDS(file = paste0("data/", cohort, "/covariates.RDS"))

reference_DNAm_corrected = correction_DNAm(bootstrap_samples     = list(rownames(covariates)),
                                           DNAm.npy              = myNorm.mdd, 
                                           i                     = 1, 
                                           pd_mdd2               = pd_mdd, 
                                           LeucocyteFraction.mdd = LeucocyteFraction.mdd, 
                                           freq                  = 0.1,
                                           CovDNAm               = CovDNAm,
                                           cohort                = cohort, 
                                           slide                 = slide)

if (is.null(reference_DNAm_corrected)){
  stop("Reference sample generated NULL output")
}

if (length(reference_DNAm_corrected) == 0){
  stop("Issue with the reference sample")
}

if (sum(colnames(reference_DNAm_corrected) == colnames(bootstrap_DNAm_corrected[[1]]))/dim(reference_DNAm_corrected)[2] !=1){
  stop("There is a different number of variables between the reference smaple and all the other bootstrap samples")
}


saveRDS(reference_DNAm_corrected , file = paste0("results/02_Preprocessing/", 
                                                 cohort, 
                                                 "/DNAm/bootstrap_DNAm_corrected_", 
                                                 experiment_name,
                                                 "_nbChunk_",
                                                 "reference",
                                                 ".RDS"))
