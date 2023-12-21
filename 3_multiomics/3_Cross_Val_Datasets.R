library(dplyr)
library(stringr)
library(DESeq2)
library(ChAMP)
library(sva)
library(SNFtool)
library(mclust)
library(caret)
library(reticulate)

################
###   mRNA   ###
################

#Covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

#Cross validation: dans Cross_val_test? ou cross_val_diff?
cv_fold = readRDS ("data/cv_fold.RDS")

#raw
all_data = readRDS ("data/all_data_mRNA.rds")

apply_deseq = function (data , covartiates) {
  
  ### vst norm 
  dds_train_adj = DESeqDataSetFromMatrix( data , 
                                          colData = covartiates[colnames(data) ,] , 
                                          design  = ~ AGE + BMI + Rin  +SEX + Polynuclear_neutrophile + Lymphocyte + GROUP)
  dds = DESeq(dds_train_adj)
  res =  results(dds ,  cooksCutoff=FALSE, independentFiltering=FALSE)
  return( as.data.frame(res) ) 
}


### filter low expressed genes
train_mRNA   = all_data [ , cv_fold[[1]]$train ]
train_mRNA_f = train_mRNA [rowMeans (train_mRNA) > 10 ,]

train_mRNA_f_female  = train_mRNA [ rownames(train_mRNA_f) , cv_fold[[1]]$train [grep ("^[0-9]" ,cv_fold[[1]]$train )] ]
train_mRNA_f_male    = train_mRNA [ rownames(train_mRNA_f) ,  cv_fold[[1]]$train [grep ("^P" ,cv_fold[[1]]$train )] ] 

diff_mRNA = function (data , folds , cov_pooled) {
  data_f = data [rowMeans (data) > 10 ,]
  diff = apply_deseq (data_f , cov_pooled)
  diffm = diff %>% filter (padj < 0.05)  %>% row.names 
  return(diffm)}

diffGenes = list()
for (i in 1: length (cv_fold) ) {
  diffGenes = append (diffGenes , list( diff_mRNA ( all_data , cv_fold [[i]] , cov_pooled) ) ) 
}

saveRDS(diffGenes , "results/diffGenesCV.rds")

#################
###   miRNA   ###
#################

raw.miRNA = read.csv2 (file = "data/miR.MDD.raw.renamed.counts.csv",
                       check.names      = FALSE,
                       stringsAsFactors = FALSE,
                       row.names        = 1)
cov.miRNA = readRDS("data/covMIRNA.mDD.RDS")

cov.miRNA = cov.miRNA[ which(cov.miRNA$ID_TGML != ""),]
raw.miRNA = raw.miRNA [ ,row.names(cov_pooled)]
raw.miRNA[is.na(raw.miRNA) ]= 0 ### NA as 0


diff_miRNA = function (data , folds , cov_pooled) {
  
  ### filter low expressed genes
  train_miRNA   = data [ , folds$train]
  Controls_train = cov_pooled[folds$train , ] %>% filter ( GROUP == "Control" ) %>% rownames
  Patient_train  = cov_pooled[folds$train , ] %>% filter ( GROUP == "MDD" ) %>% rownames

  # filter miRs : a mir is retained only if present with more then 1
  # reads in more then 60% of either patients or controls 
  
  x_train = rowSums(train_miRNA[,Patient_train]  > 1) >= 0.6*length(Patient_train)
  y_train = rowSums(train_miRNA[,Controls_train] > 1) >= 0.6*length(Controls_train)
  raw.miRNA.f.train = train_miRNA [x_train | y_train , ]
  
  
  ### vst norm 
  dds_train_adj = DESeqDataSetFromMatrix( raw.miRNA.f.train , 
                                          colData = cov_pooled[colnames(raw.miRNA.f.train) ,] , 
                                          design  = ~ batchmiRNA + AGE + BMI + Rin.miR  + SEX + Polynuclear_neutrophile + Lymphocyte + GROUP)
  dds = DESeq(dds_train_adj)
  res =  results(dds ,  cooksCutoff=FALSE, independentFiltering=FALSE)
  
  diffMIRs = res %>% as.data.frame %>% filter (padj < 0.05) %>% rownames
  return( diffMIRs  ) 
  
  
}

diffmiRS = list()
for (i in 1: length (cv_fold) ) {
  diffmiRS = append (diffmiRS , list( diff_miRNA ( raw.miRNA , cv_fold [[i]] , cov_pooled) ) ) 
}
saveRDS(diffmiRS , "results/diffmiRSCV.rds")

################
###   DNAm   ###
################

pd_mdd_bin = readRDS("D:/CC/Scripts/R/gitHub_MDD/multi_omics/data/pd_mdd_bin.RDS") # pd file containes metadata of samples
myNorm.mdd = readRDS("D:/CC/Scripts/R/gitHub_MDD/multi_omics/data/myNorm.mdd.RDS") # normalised beta-values of probes
LeucocyteFraction.mdd = readRDS("D:/CC/Scripts/R/gitHub_MDD/data/LeucocyteFraction.mdd.rds")
td = readRDS ("D:/CC/Scripts/R/gitHub_MDD/multi_omics/data/pd_mdd_female.RDS")
female_mdd = readRDS("D:/CC/Scripts/R/gitHub_MDD/multi_omics/data/cov.icm.RDS")
male_mdd = readRDS("D:/CC/Scripts/R/gitHub_MDD/multi_omics/data/cov.igbmc.RDS")

name_grp = rbind(male_mdd  [,c("Group","Name")] , female_mdd  [,c("Group","Name")]  )
rownames(name_grp) = rownames(name_grp) %>% gsub( "_L[0-9]+", "" , .)

pd_mdd = pd_mdd_bin [match(name_grp [ rownames(cov_pooled) , "Name"] , pd_mdd_bin$Sample_ID ),]
pd_mdd$Name = rownames(cov_pooled)
rownames(pd_mdd) =  pd_mdd$Name
pd_mdd2 = pd_mdd
pd_mdd2$Sex = pd_mdd2$Sex %>% as.factor %>% as.double
pd_mdd2$Array = pd_mdd2$Array %>% as.factor %>% as.double
pd_mdd2$BMI.bin = pd_mdd2$BMI.bin %>% as.factor %>% as.double
pd_mdd2$Slide = pd_mdd2$Slide %>% as.factor %>% as.double

DNAm.npy = myNorm.mdd
colnames(DNAm.npy) = pd_mdd [match ( colnames(DNAm.npy) , pd_mdd$Sample_Name  ) ,] %>% rownames
DNAm.npy = DNAm.npy [ , rownames(cov_pooled) ]
DNAm.npy = t(DNAm.npy)

train_nms = paste0( "CvDiff_" , rep (1:5 , each = 5) ,"_p" , 1:5 ,"_train_DNA_m.RDS" )
test_nms = paste0( "CvDiff_" , rep (1:5 , each = 5) ,"_p" , 1:5 ,"_test_DNA_m.RDS" )
DNAmDiffCv_nms = paste0( "DiffDNAmCV_" , rep (1:5 , each = 5) ,"_p" , 1:5 ,"_train_DNA_m.RDS" )

DNAmCVS = lapply ( 1:length(train_nms) , function (x) readRDS ( paste0("results/CV/" , train_nms [[x]]) ) )
corrected_m.value.mdd.bin.nopreservation = lumi::beta2m(DNAmCVS[[1]])

