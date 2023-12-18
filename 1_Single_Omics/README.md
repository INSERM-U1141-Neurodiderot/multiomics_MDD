---
title: "MDD Single Modality analyse"
author: "U1141"
date: "2023-12-07"
output: html_document
---
# MDD_Single_Modality
Case-Control single omics analysis for DNAm, mRNA and miRNA.

# DNAm

```{r cars}
library(sva)
library(ChAMP)
library(lumi)
```

#### Data
```{r cars}
# data = readRDS('data_DNAm.rds')
beta = data$beta
pd_mdd_bin = data$pd_mdd_bin
sample_sheet_covariates = data$sample_sheet_covariates
sample_sheet = data$sample_sheet
rm(data)
```
Raw data used
```{r cars}
dim(beta)
```
Covariates
```{r cars}
summary(pd_mdd_bin)
```
```{r cars}
summary(sample_sheet_covariates)
```
Technical data
```{r cars}
summary(sample_sheet)
```

Cell count (for data correction)
```{r cars}
LeucocyteFraction.mdd = readRDS("LeucocyteFraction.mdd.rds")
```

Quality control
```{r cars}
CpG.GUI(CpG = rownames(beta), arraytype = "EPIC")
```
```{r cars}
champ.QC(beta = beta, pheno = sample_sheet$Sample_Group,
         resultsDir = "QC_beforeNormalisation", Feature.sel = "SVD")
```


#### Normalization
Data normalisation using ChAMP
```{r cars}
myNorm = champ.norm(beta = beta, method = "BMIQ", plotBMIQ = T,
                    arraytype = "EPIC", resultsDir = "Normalisation", cores = 2)
```
QUality after normalization
```{r cars}
champ.QC(beta = as.matrix(myNorm), pheno = sample_sheet$Sample_Group, resultsDir = "QC_afterNorm")
```

#### Data Correction
Creating model
```{r cars}
mod = model.matrix(~1, data = pd_mdd_bin)
```
Correction
```{r cars}
bat = myNorm[, pd_mdd_bin$Sample_Name]
for (col in c('Slide', 'Array', 'Sex', 'Age_bin', 'BMI.bin')){
  print(col)
  batch = pd_mdd_bin[, col]
  i = 0
  for (n in names(table(batch))){
    batch = ifelse(batch == n, i, batch)
    i = i + 1
  }
  bat = ComBat(dat = bat, batch = batch, mod = mod)
}
```
#### Cell count correction
Correction
```{r cars}
blood = LeucocyteFraction.mdd[colnames(bat),]
beta.lm = apply(bat, 1, function(x) {
  LeucocyteFraction.mdd[colnames(bat),] = blood 
  lm(x~CD4+CD8+MO+B+NK+GR, data = blood)
})
```
Residuals
```{r cars}
residuals = t(sapply(beta.lm, function(x){residuals(summary(x))}))
colnames(residuals) = colnames(bat)
plot(residuals )
```
Correction
```{r cars}
b.value.mdd.bin = residuals + matrix(apply(bat, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
b.value.mdd.bin[b.value.mdd.bin >= 1] = 0.99999999
b.value.mdd.bin[b.value.mdd.bin <= 0] = 0.00000001

corrected_m.value.mdd.bin.nopreservation = lumi::beta2m(b.value.mdd.bin)
```
#### DMP analysis
Detecting differentially methylated probes on Female cohort
```{r cars}
female = corrected_m.value.mdd.bin.nopreservation[, as.character(unlist(sample_sheet_covariates[sample_sheet_covariates$Sex == 'F', 'Sample_Name']))]
dmp_female = champ.DMP(beta = female, pheno = pd_mdd_bin[pd_mdd_bin$Sex == 'F', 'Sample_Group'],
                       compare.group = c("control", "mdd"), adjPVal = 1,
                       adjust.method = "BH", arraytype = "EPIC")
dmp_female = as.data.frame(dmp_female$control_to_mdd)
```
Male cohort
```{r cars}
male = corrected_m.value.mdd.bin.nopreservation[, as.character(unlist(sample_sheet_covariates[sample_sheet_covariates$Sex == 'M', 'Sample_Name']))]
dmp_male = champ.DMP(beta = male, pheno = pd_mdd_bin[pd_mdd_bin$Sex == 'M', 'Sample_Group'],
                     compare.group = c("control", "mdd"), adjPVal = 1,
                     adjust.method = "BH", arraytype = "EPIC")
dmp_male = as.data.frame(dmp_male$control_to_mdd)
```
Significant DMP (p-value < 0.05)
```{r cars}
dmp_female_sig = dmp_female[dmp_female$P.Value < 0.05,]
dmp_male_sig = dmp_male[dmp_male$P.Value < 0.05,]
dim(dmp_female_sig)
dim(dmp_male_sig)
```

# mRNA
```{r cars}
library(variancePartition)
library(DESeq2)
library(dplyr)
```
#### Data
```{r cars}
data = readRDS('data.mRNA.rds')
```
```{r cars}
raw.data.Female = data$raw.data.Female
raw.data.Male = data$raw.data.Male
cov.Female = data$cov.Female
cov.Male = data$cov.Male
rm(data)
```
#### Quality control
Minimum read counts for Male
```{r cars}
dim(cov.Male)
Male_filter = raw.data.Male[which(rowMeans(raw.data.Male) > 10), ]
cov.Male = cov.Male[rownames(cov.Male) %in% colnames(Male_filter),]
dim(cov.Male)
```
Female
```{r cars}
dim(cov.Female)
Female_filter = raw.data.Female[which(rowMeans(raw.data.Female) > 10), ]
cov.Female = cov.Female[rownames(cov.Female) %in% colnames(Female_filter),]
dim(cov.Female)
```
#### Analysis
Analysis of differential expression
```{r cars}
Analys_diff = function (data, covariates) {
  dds = DESeqDataSetFromMatrix(countData = data,
                               colData = covariates[colnames(data),],
                               design = ~ Age + BMI + RIN + Polynuclear_neutrophile + Lymphocyte + Group)
  dds = DESeq(dds)
  
  res =  results(dds, cooksCutoff = F, independentFiltering = F)
  cat(summary(res))
  cat('\n')
  
  res = as.data.frame(res)
  res = res[order(res$log2FoldChange),]
  
  return(res)
}

analys_diff_Male = Analys_diff(data = Male_filter, covariates = cov.Male)
analys_diff_Female = Analys_diff(Female_filter, cov.Female)
```
Significant differentialy expressed genes (DEG)
```{r cars}
DEG_male = analys_diff_Male[analys_diff_Male$padj < 0.05,]
DEG_female = analys_diff_Female[analys_diff_Female$padj < 0.05,]
dim(DEG_male)
dim(DEG_female)
```

# miRNA
```{r cars}
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
```

#### Data
Covariables
```{r cars}
cov.miRNA = readRDS(file = "data/covMIRNA.mDD.RDS")
cov.miRNA = cov.miRNA[ which(cov.miRNA$ID_TGML != ""),]
liste_col = c("Age", "BMI", "Rin.miR", "Rin.mRNA",
              "Polynuclear_neutrophile", "Leukocyte",
              "Lymphocyte", "Platelet")
cov.miRNA[,liste_col] = apply(cov.miRNA[, liste_col], 2, as.numeric)
cov.miRNA$Sex.num = as.numeric(as.factor(cov.miRNA[, c("Sex")]))

cov_male = cov.miRNA[cov.miRNA$Sex == 'M',]
cov_female =cov.miRNA[cov.miRNA$Sex == 'F',]
```
Raw data
```{r cars}
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
  
  #Filter : 60% reads to 1 min for case or ctr
  raw = raw [rowSums(raw[,Patient]  > 1) >= 0.6*length(Patient) |
               rowSums(raw[,Controls] > 1) >= 0.6*length(Controls), ]
  return(raw)
}
raw.miRNA.male = format_mir(raw.miRNA[, rownames(cov.miRNA[cov.miRNA$Sex == 'M',])], cov_male)
raw.miRNA.female = format_mir(raw.miRNA[, rownames(cov.miRNA[cov.miRNA$Sex == 'F',])], cov_female)
```
#### Imputation
Imputing missing values
```{r cars}
imputer = function(data, groupe, col){
  data[which(data$Group == groupe & is.na(data[, col])) , col ] =
    median(na.omit(data[which(data$Group == groupe & !is.na(data[, col])), col]))
  return(data)
}
table(is.na(cov_male))
table(is.na(cov_female))
for (groupe in c('Patient', 'Control')){
  for (col in liste_col){
    cov_male = imputer(cov_male, groupe, col)
    cov_female = imputer(cov_female, groupe, col)
  }
}
```
#### Analysis
```{r cars}
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
```
Significant miRNA
```{r cars}
mir_male = res_male$res[res_male$res $padj < 0.05,]
mir_female = res_female$res[res_female$res$padj < 0.05,]
dim(mir_male)
dim(mir_female)
```







