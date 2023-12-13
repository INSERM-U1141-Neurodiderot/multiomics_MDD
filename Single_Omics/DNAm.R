library(sva)
library(ChAMP)

####################
###     DATA     ###
####################

data = readRDS('data/data_DNAm.rds')
beta = data$beta
pd_mdd_bin = data$pd_mdd_bin
sample_sheet_covariates = data$sample_sheet_covariates
sample_sheet = data$sample_sheet
rm(data)

#Cell count
LeucocyteFraction.mdd = readRDS("data/LeucocyteFraction.mdd.rds")

#Quality
CpG.GUI(CpG = rownames(beta), arraytype = "EPIC")
champ.QC(beta = beta, pheno = sample_sheet$Sample_Group,
         resultsDir = "QC_beforeNormalisation", Feature.sel = "SVD")

#############################
###     Normalisation     ###
#############################
myNorm = champ.norm(beta = beta, method = "BMIQ", plotBMIQ = T,
                    arraytype = "EPIC", resultsDir = "Normalisation", cores = 2)
champ.QC(beta = as.matrix(myNorm), pheno = sample_sheet$Sample_Group, resultsDir = "QC_afterNorm")

###########################
###   Data Correction   ###
###########################

mod = model.matrix(~1, data = pd_mdd_bin)
bat = myNorm[, pd_mdd_bin$Sample_Name]
for (col in c('Slide', 'Array', 'Sex', 'Age_bin', 'BMI.bin')){
  batch = pd_mdd_bin[, col]
  i = 0
  for (n in names(table(batch))){
    batch = ifelse(batch == n, i, batch)
    i = i + 1
  }
  bat = ComBat(dat = bat, batch = batch, mod = mod)
  print(col)
}

#################################
###   Cell count correction   ###
#################################
blood = LeucocyteFraction.mdd[colnames(bat),]
beta.lm = apply(bat, 1, function(x) {
  LeucocyteFraction.mdd[colnames(bat),] = blood 
  lm(x~CD4+CD8+MO+B+NK+GR, data = blood)
})

residuals = t(sapply(beta.lm, function(x){residuals(summary(x))}))
colnames(residuals) = colnames(bat)

b.value.mdd.bin = residuals + matrix(apply(bat, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
b.value.mdd.bin[b.value.mdd.bin >= 1] = 0.99999999
b.value.mdd.bin[b.value.mdd.bin <= 0] = 0.00000001

corrected_m.value.mdd.bin.nopreservation = lumi::beta2m(b.value.mdd.bin)

########################
###   DMP analysis   ###
########################
female = corrected_m.value.mdd.bin.nopreservation[, as.character(unlist(sample_sheet_covariates[sample_sheet_covariates$Sex == 'F', 'Sample_Name']))]
male = corrected_m.value.mdd.bin.nopreservation[, as.character(unlist(sample_sheet_covariates[sample_sheet_covariates$Sex == 'M', 'Sample_Name']))]
dmp_male = champ.DMP(beta = male, pheno = pd_mdd_bin[pd_mdd_bin$Sex == 'M', 'Sample_Group'],
                     compare.group = c("control", "mdd"), adjPVal = 1,
                     adjust.method = "BH", arraytype = "EPIC")
dmp_female = champ.DMP(beta = female, pheno = pd_mdd_bin[pd_mdd_bin$Sex == 'F', 'Sample_Group'],
                       compare.group = c("control", "mdd"), adjPVal = 1,
                       adjust.method = "BH", arraytype = "EPIC")
dmp_male = as.data.frame(dmp_male$control_to_mdd)
dmp_female = as.data.frame(dmp_female$control_to_mdd)
dmp_male_sig = dmp_male[dmp_male$adj.P.Val < 0.05,]
dmp_female_sig = dmp_female[dmp_female$adj.P.Val < 0.05,]

saveRDS(list(dmp_male = dmp_male, dmp_female = dmp_female), 'res_single_DNAm.rds')
