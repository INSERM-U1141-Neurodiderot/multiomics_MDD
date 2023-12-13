library(variancePartition)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(multiMiR)
library(ggpubr)
library(FactoMineR)
library(ggcorrplot)
library(dendextend)
library(factoextra)
library(corrr)
library(corrplot)

################
###   Data   ###
################

cov.miRNA = readRDS(file = "data/covMIRNA.mDD.RDS")
cov.miRNA = cov.miRNA[ which(cov.miRNA$ID_TGML != ""),]

liste_col = c("Age", "BMI", "Rin.miR", "Rin.mRNA",
              "Polynuclear_neutrophile", "Leukocyte",
              "Lymphocyte", "Platelet")
cov.miRNA[,liste_col] = apply(cov.miRNA[, liste_col], 2, as.numeric)
cov.miRNA$Sex.num = as.numeric(as.factor(cov.miRNA[, c("Sex")]))

raw.miRNA = read.csv2(file = "data/miR.MDD.raw.renamed.counts.csv",
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = 1)

format_mir = function(raw, cov){
  #NA à 0
  raw[is.na(raw)] = 0
  
  #MDD / CTR
  Controls = rownames(cov[cov$Group == "Control",])
  Patient = rownames(cov[cov$Group == "Patient",])
  
  #Filtre : 60% de reads à 1 minium pour patient ou ctr
  raw = raw [rowSums(raw[,Patient]  > 1) >= 0.6*length(Patient) |
               rowSums(raw[,Controls] > 1) >= 0.6*length(Controls), ]
  return(raw)
}
cov_male = cov.miRNA[cov.miRNA$Sex == 'M',]
cov_female =cov.miRNA[cov.miRNA$Sex == 'F',]
raw.miRNA.male = format_mir(raw.miRNA[, rownames(cov.miRNA[cov.miRNA$Sex == 'M',])], cov_male)
raw.miRNA.female = format_mir(raw.miRNA[, rownames(cov.miRNA[cov.miRNA$Sex == 'F',])], cov_female)

####################
###   Imputing   ###
####################

imputer = function(data, groupe, col){
  data[which(data$Group == groupe & is.na(data[, col])) , col ] =
    median(na.omit(data[which(data$Group == groupe & !is.na(data[, col])), col]))
  return(data)
}
for (groupe in c('Patient', 'Control')){
  for (col in liste_col){
    cov_male = imputer(cov_male, groupe, col)
    cov_female = imputer(cov_female, groupe, col)
  }
}


####################
###   Analysis   ###
####################

analyse_mir = function(data, cov){
  dds.Mir = DESeqDataSetFromMatrix(countData = data[, rownames(cov)] ,
                                   colData = cov,
                                   design = ~ Rin.miR + batchmiRNA + Group)
  res = DESeq(dds.Mir)
  res =  results(res, cooksCutoff = F, independentFiltering = F)
  res = as.data.frame(res)
  res = res[order(res$pvalue),]
  
  # Normalized
  dds.Mir = estimateSizeFactors(dds.Mir)
  dds.Mir.vst = varianceStabilizingTransformation(dds.Mir)
  
  dds.Mir = estimateDispersions(dds.Mir)
  dds.Mir = nbinomWaldTest(dds.Mir)
  res.Mir = results(dds.Mir, cooksCutoff = T)
  
  miR.corrected.Data.vst = t(lm(t(assay(dds.Mir.vst))~ Rin.miR + BMI + Age + batchmiRNA,
                                data = cov[colnames(assay(dds.Mir.vst)),])$residuals)
  miR.corrected.Data.vst = round((miR.corrected.Data.vst + abs(min(miR.corrected.Data.vst)))*10)
  dds.Mir = DESeqDataSetFromMatrix(countData = miR.corrected.Data.vst[, rownames(cov)] ,
                                   colData = cov,
                                   design = ~ Rin.miR + batchmiRNA + Group)
  res.norm = DESeq(dds.Mir)
  res.norm =  results(res.norm, cooksCutoff = F, independentFiltering = F)
  res.norm = as.data.frame(res.norm)
  res.norm = res.norm[order(res.norm$pvalue),]
  
  return(list(res = res, res.norm = res.norm))
}
res_female = analyse_mir(raw.miRNA.female, cov_female)
res_male = analyse_mir(raw.miRNA.male, cov_male)

saveRDS(list(miR_male = res_male, miR_female = res_female), 'res_single_miRNA.rds')
