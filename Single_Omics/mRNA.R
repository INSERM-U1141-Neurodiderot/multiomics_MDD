library(variancePartition)
library(DESeq2)
library(dplyr)
library(ggplot2)

################
###   Data   ###
################
data = readRDS('data/data.mRNA.rds')
raw.data.Female = data$raw.data.Female
raw.data.Male = data$raw.data.Male
cov.Female = data$cov.Female
cov.Male = data$cov.Male

rm(data)

###################
###   Quality   ###
###################

Male_CTRL = (row.names(cov.Male[cov.Male$Group == "Control",]))
Male_MDD = (row.names(cov.Male[cov.Male$Group == "Patient",]))

Female_CTRL = (row.names(cov.Female[cov.Female$Group == "Control",]))
Female_MDD = (row.names(cov.Female[cov.Female$Group == "Patient",]))

Male_filter = raw.data.Male[which(rowMeans(raw.data.Male) > 10), ]
Female_filter = raw.data.Female[which(rowMeans(raw.data.Female) > 10), ]

cov.Male = cov.Male[rownames(cov.Male) %in% colnames(Male_filter),]
cov.Female = cov.Female[rownames(cov.Female) %in% colnames(Female_filter),]

####################
###   Analysis   ###
####################
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


saveRDS(list(diff_male = analys_diff_Male, diff_female = analys_diff_Female), 'res_single_mRNA.rds')
