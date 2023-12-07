---
title: "MDD_codes"
author: "U1141"
date: "2023-12-07"
output: html_document
---
---
title: "MDD_Single_Modality"
author: "U1141"
date: "2023-12-06"
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

# Rank-rank hypergeometric overlap (RRHO) analysis
```{r cars}
library(dplyr)
require(parallel)
require(VennDiagram)
require(RRHO)
```

For each omic, results from differential analysis in each sex were ranked based on the following metric: -log10(p-value) x sign(log2 Fold Change). Then, the RRHO2 function was applied to the 2 lists at default parameters (with step size equal to the square root of the list length)

#### Data
Using data generated by DEG analysis
```{r}
data = readRDS('res_single_mRNA.rds')

female = data$diff_female
male = data$diff_male
rm(data)
```

DDE: -log10(P-value)xFoldChangeSign
```{r }
female$DDE = -log10(female$pvalue) * ifelse(female$log2FoldChange > 0, 1, -1)
female$Genes = rownames(female)
female = female[, c('Genes', 'DDE')]

male$DDE = -log10(male$pvalue) * ifelse(male$log2FoldChange > 0, 1, -1)
male$Genes = rownames(male)
male = male[, c('Genes', 'DDE')]

male = male[male$Genes %in% female$Genes,]
female = female[female$Genes %in% male$Genes,]
```

## Initialization
```{r }
stepsize = RRHO:::defaultStepSize(male, female)
male <- male[order(male[, 2], decreasing = TRUE), ]
female <- female[order(female[, 2], decreasing = TRUE), ]
nmale <- length(male[, 1])
nfemale <- length(female[, 1])
N <- max(nmale, nfemale)

labels = c("male", "female")
boundary = 0.1
```

Computing the overlaps between two *numeric* lists:
```{r }
numericListOverlap<- function(sample1, sample2, stepsize, method="hyper", tol = 0.5, offset = 1 , mcores = ncores){
  n<- length(sample1)
  
  overlap_hyper <- function(a,b) {
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))    
    signs<- 1L
    log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)    
    #log.pval[is.na(log.pval)]<-0
    
    return(c(counts=count, 
             log.pval=as.numeric(log.pval),
             signs=as.integer(signs)
    ))    
  }
  
  overlap_fisher <- function(a,b) {
    s1 <- sample1[1:a]
    s2 <- sample2[1:b]
    lenA <- as.integer(sum(as.numeric(s1 %in% s2))) 
    lenB <- length(s1)
    lenC <- length(s2)
    Odds<-((lenA+offset)*(n-lenB-lenC+lenA+offset))/((lenC-lenA+offset)*(lenB-lenA+offset))
    logOdds <- log(abs(Odds))*sign(Odds)
    signs<- 1L
    
    return(c(counts=lenA, 
             log.pval=as.numeric(logOdds),
             signs=as.integer(signs)
    ))    
  }
  
  indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize))
  if(method=="hyper"){
    #overlaps<- apply(indexes, 1, function(x) overlap_hyper(x['i'], x['j']))
    overlaps<- parallel::mcmapply(overlap_hyper , a = indexes[ ,1] , b = indexes[ ,2] , mc.cores = mcores)
  } else if(method=="fisher"){
    #overlaps<- apply(indexes, 1, function(x) overlap_fisher(x['i'], x['j']))
    overlaps<- parallel::mcmapply(overlap_fisher , a = indexes[ ,1] , b = indexes[ ,2] , mc.cores = mcores)
  }
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows)  
  matrix.signs<- matrix(overlaps['signs',], ncol=nrows)  
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals, signs = matrix.signs))
}
```

Boundaries
```{r }
.hypermat_normal<- numericListOverlap(male[, 1], female[, 1], stepsize, method='hyper' , mcores = 1)
hypermat_normal<- .hypermat_normal$log.pval

.hypermat_flipX <- numericListOverlap(rev(male[, 1]), female[, 1], stepsize, method='hyper', mcores = 1)
hypermat_flipX <- .hypermat_flipX$log.pval

stepmale <- seq(1, nmale, stepsize)
stepfemale <- seq(1, nfemale, stepsize)

len_male <- length(stepmale)
len_female <- length(stepfemale)

lenStrip_m <- round(len_male*boundary)
lenStrip_f <- round(len_female*boundary)

boundary_m <- sum(male[stepmale,2] > 0)
boundary_f <- sum(female[stepfemale,2] > 0)
```
Initialize matrix
```{r }
hypermat <- matrix(NA, nrow = nrow(hypermat_normal) + lenStrip_m,
                   ncol = ncol(hypermat_normal) + lenStrip_f)
# d1d2, quadrant I
hypermat[lenStrip_m + (boundary_m+1):len_male, lenStrip_f + (boundary_f+1):len_female] <- hypermat_normal[(boundary_m+1):len_male, (boundary_f+1):len_female]

# u1d2, quadrant II
hypermat[1:boundary_m, lenStrip_f + (boundary_f+1):len_female] <- hypermat_flipX[len_male:(len_male - boundary_m + 1),(boundary_f+1):len_female]
  
# u1u2, quadrant III
hypermat[1:boundary_m, 1:boundary_f] <- hypermat_normal[1:boundary_m,1:boundary_f]
  
# u1d2, quadrant IV
hypermat[lenStrip_m + (boundary_m+1):len_male, 1:boundary_f] <- hypermat_flipX[(len_male - boundary_m):1,1:boundary_f]

hypermat <- hypermat * log10(exp(1))
```

dd: down in 1 and down in 2
```{r }
maxind.dd <- which(max(hypermat[lenStrip_m + (boundary_m+1):len_male , lenStrip_f + (boundary_f+1):len_female],
                       na.rm = TRUE) == hypermat, arr.ind = TRUE)
maxind.dd <- maxind.dd[maxind.dd[,1]>=lenStrip_m + (boundary_m+1) & maxind.dd[,1]<=lenStrip_m +len_male & 
                         maxind.dd[,2]>=lenStrip_f + (boundary_f+1) & maxind.dd[,2]<=lenStrip_f + len_female,]

if(!is.null(dim(maxind.dd))){
  maxind.dd <- maxind.dd[1, ]
}

indmale.dd <- seq(1, nmale, stepsize)[maxind.dd[1] - lenStrip_m]
indfemale.dd <- seq(1, nfemale, stepsize)[maxind.dd[2] - lenStrip_f]
gene_male_dd <- male[indmale.dd:nmale, 1]
gene_female_dd <- female[indfemale.dd:nfemale, 1]
gene_list_overlap_dd <- intersect(gene_male_dd,
                                  gene_female_dd)
genelist_dd <- list(gene_male_dd=gene_male_dd, 
                    gene_female_dd=gene_female_dd,
                    gene_list_overlap_dd=gene_list_overlap_dd
)
```
uu: up in 1 and up in 2
```{r }
maxind.uu <- which(max(hypermat[1:boundary_m, 1:boundary_f],
                       na.rm = TRUE) == hypermat, arr.ind = TRUE)
maxind.uu <- maxind.uu[maxind.uu[,1]>=1 & maxind.uu[,1]<=boundary_m & maxind.uu[,2]>=1 & maxind.uu[,2]<=boundary_f,]
if(!is.null(dim(maxind.uu))){
  maxind.uu <- maxind.uu[1, ]
}

indmale.uu <- seq(1, nmale, stepsize)[maxind.uu[1]]
indfemale.uu <- seq(1, nfemale, stepsize)[maxind.uu[2]]
gene_male_uu <- male[1:indmale.uu, 1]
gene_female_uu <- female[1:indfemale.uu, 1]
gene_list_overlap_uu <- intersect(gene_male_uu,
                                  gene_female_uu)
genelist_uu <- list(gene_male_uu=gene_male_uu, 
                    gene_female_uu=gene_female_uu,
                    gene_list_overlap_uu=gene_list_overlap_uu
)
```

ud: up in 1 and down in 2
```{r }
maxind.ud <- which(max(hypermat[1:boundary_m, lenStrip_f + (boundary_f+1):len_female],
                       na.rm = TRUE) == hypermat, arr.ind = TRUE)
#
maxind.ud <- maxind.ud[maxind.ud[,1]>=1 & maxind.ud[,1]<=boundary_m & maxind.ud[,2]>= lenStrip_f + (boundary_f+1) & maxind.ud[,2]<=lenStrip_f + len_female,]
if(!is.null(dim(maxind.ud))){
  maxind.ud <- maxind.ud[1, ]
}

indmale.ud <- seq(1, nmale, stepsize)[maxind.ud[1]]
indfemale.ud <- seq(1, nfemale, stepsize)[maxind.ud[2] - lenStrip_f]
gene_male_ud <- male[1:indmale.ud, 1]
gene_female_ud <- female[indfemale.ud:nfemale, 1]
gene_list_overlap_ud <- intersect(gene_male_ud,
                                  gene_female_ud)
genelist_ud <- list(gene_male_ud=gene_male_ud, 
                    gene_female_ud=gene_female_ud,
                    gene_list_overlap_ud=gene_list_overlap_ud
)
```

du: down in 1 and up in 2
```{r }
maxind.du <- which(max(hypermat[lenStrip_m + (boundary_m+1):len_male, 1:boundary_f],
					 na.rm = TRUE) == hypermat, arr.ind = TRUE)
#
maxind.du <- maxind.du[maxind.du[,1]>=lenStrip_m + (boundary_m+1) & maxind.du[,1]<=lenStrip_m + len_male & maxind.du[,2]>=1 & maxind.du[,2]<=boundary_f,]
if(!is.null(dim(maxind.du))){
maxind.du <- maxind.du[1, ]
}

indmale.du <- seq(1, nmale, stepsize)[maxind.du[1] - lenStrip_m]
indfemale.du <- seq(1, nfemale, stepsize)[maxind.du[2]]
gene_male_du <- male[indmale.du:nmale, 1]
gene_female_du <- female[1:indfemale.du, 1]
gene_list_overlap_du <- intersect(gene_male_du,
								gene_female_du)
genelist_du <- list(gene_male_du=gene_male_du, 
				  gene_female_du=gene_female_du,
				  gene_list_overlap_du=gene_list_overlap_du
)
```

## Best values
get best log10(p-val) of each quadrent
```{r }
max_na_row = which (is.na(hypermat [1,])) %>% max 
max_na_col = which (is.na(hypermat [,1])) %>% max 
min_na_row = which (is.na(hypermat [1,])) %>% min 
min_na_col = which (is.na(hypermat [,1])) %>% min 


best_UU = (as.numeric(unlist(hypermat [1:min_na_row, 1:min_na_col] )) %>% na.omit %>% max )
best_DD = (as.numeric(unlist(hypermat [max_na_row:nrow(hypermat) ,max_na_col:ncol(hypermat)] )) %>% na.omit %>% max )
best_UD = (as.numeric(unlist(hypermat [max_na_row:nrow(hypermat), 1:min_na_col] )) %>% na.omit %>% max ) 
best_DU = (as.numeric(unlist(hypermat [1:min_na_row , max_na_col:ncol(hypermat)] )) %>% na.omit %>% max )

best_UU
best_DD
best_UD
best_DU
```
## Heatmap

```{r }
maximum <- max(hypermat,na.rm=TRUE)
minimum <- min(hypermat,na.rm=TRUE)

color.bar <- function(lut, min, max=-min, 
                      nticks=11, 
                      ticks=seq(min, max, len=nticks), 
                      title='') {
  scale  <- (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', 
       xaxt='n', xlab='', yaxt='n', ylab='')
  mtext(title,2,2.3, cex=0.8)
  axis(2, round(ticks,0), las=1,cex.lab=0.8)
  for (i in 1:(length(lut)-1)) {
    y  <- (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
  

jet.colors  <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colorGradient <- jet.colors(101)

layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  
breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
image(hypermat, col = colorGradient,breaks=breaks,vaxes = FALSE)
  
mtext('female',2,0.5)
mtext('male',1,0.5)
  
color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = '-log10(P-value)')
invisible(hypermat)

```






