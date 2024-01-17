library(dplyr)
library(stringr)
library(DESeq2)
library(ChAMP)
library(sva)

### input: covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

### input: cv folds
cv_fold = readRDS("results/1_CrossValidation/cv_fold.RDS")
cv_fold_female = readRDS("results/1_CrossValidation/cv_fold_female.RDS")
cv_fold_male = readRDS("results/1_CrossValidation/cv_fold_male.RDS")

################
###   mRNA   ###
################

data = readRDS('data/data_mRNA.rds')

correct_mRNA = function (data , folds , cov_pooled, All = TRUE) {
        
    ### filter low expressed genes
    train_mRNA   = data [ , folds$train]
    train_mRNA_f = train_mRNA [rowMeans (train_mRNA) > 10 ,]

    test_mRNA_f  = data [ rownames(train_mRNA_f) , folds$test]
    
    
    ### vst norm 
    dds_train_adj = DESeqDataSetFromMatrix( (train_mRNA_f) , colData=cov_pooled[folds$train,] , design=~1)
    vst_train_adj = vst(dds_train_adj)

    dds_test_adj = DESeqDataSetFromMatrix( (test_mRNA_f) , colData=cov_pooled[folds$test,] , design=~1)
    vst_test_adj = vst(dds_test_adj)
    
    ### Correct data 
    if (All == TRUE){
        covmRNA  = c("BMI",  "Rin", "AGE", "SEX" ,"Polynuclear_neutrophile", "Lymphocyte")
        
        corrected_mRNA_train_model = (lm(t(assay(vst_train_adj))~ Rin + BMI + AGE + SEX + Polynuclear_neutrophile + Lymphocyte, 
                           data=cov_pooled[colnames(assay(vst_train_adj)) , covmRNA]) )
        
        corrected_mRNA_train = corrected_mRNA_train_model$residuals
        
        corrected_mRNA_test_p = predict (corrected_mRNA_train_model ,  cbind( t(assay(vst_test_adj)) %>% as.data.frame , cov_pooled[colnames(assay(vst_test_adj)) , covmRNA] ) )
    
        corrected_mRNA_test  =  t(assay(vst_test_adj)) %>% as.data.frame - corrected_mRNA_test_p
    }else{
        covmRNA  = c("BMI",  "Rin", "AGE" ,"Polynuclear_neutrophile", "Lymphocyte")
        
        corrected_mRNA_train_model = (lm(t(assay(vst_train_adj))~ Rin + BMI + AGE  + Polynuclear_neutrophile + Lymphocyte, 
                           data=cov_pooled[colnames(assay(vst_train_adj)) , covmRNA]) )
    
        corrected_mRNA_train = corrected_mRNA_train_model$residuals
        
        corrected_mRNA_test_p = predict (corrected_mRNA_train_model ,  cbind( t(assay(vst_test_adj)) %>% as.data.frame , cov_pooled[colnames(assay(vst_test_adj)) , covmRNA] ) )
    
        corrected_mRNA_test  =  t(assay(vst_test_adj)) %>% as.data.frame - corrected_mRNA_test_p 
    }
        
    

    
    return ( list ( corrected_mRNA_train = corrected_mRNA_train ,
                    corrected_mRNA_test  = corrected_mRNA_test) )
}

cv_mRNA_corr = correct_mRNA (data , cv_fold[[1]] , cov_pooled  )
cv_mRNA_f_corr = correct_mRNA (data , cv_fold_female[[1]] , cov_pooled , All = F)
cv_mRNA_m_corr = correct_mRNA (data , cv_fold_male[[1]] , cov_pooled , All = F)

#For correcting all folds:
# cv_mRNA_corr = lapply (cv_fold , function (x) correct_mRNA (data , x , cov_pooled  )  )
# cv_mRNA_f_corr = lapply (cv_fold_female , function (x) correct_mRNA (data , x , cov_pooled , All = FALSE  )  )
# cv_mRNA_m_corr = lapply (cv_fold_male , function (x) correct_mRNA (data , x , cov_pooled , All = FALSE  )  )

saveRDS (cv_mRNA_corr , file = "results/2_PreProcessing/cv_mRNA_corr.RDS")
saveRDS (cv_mRNA_f_corr , file = "results/2_PreProcessing/cv_mRNA_f_corr.RDS")
saveRDS (cv_mRNA_m_corr , file = "results/2_PreProcessing/cv_mRNA_m_corr.RDS")

#################
###   miRNA   ###
#################

raw.miRNA = readRDS('data/data_miRNA.rds')

correct_miRNA = function (data , folds , cov_pooled, All = T) {
        
    ### filter low expressed genes
    train_miRNA   = data [ , folds$train]
    
    Controls_train = cov_pooled[folds$train , ] %>% filter ( GROUP == "Control" ) %>% rownames
    Patient_train  = cov_pooled[folds$train , ] %>% filter ( GROUP == "MDD" ) %>% rownames
    
    # filter miRs : a mir is retained only if present with more then 1
    # reads in more then 60% of either patients or controls 

    x_train = rowSums(train_miRNA[,Patient_train]  > 1) >= 0.6*length(Patient_train)
    y_train = rowSums(train_miRNA[,Controls_train] > 1) >= 0.6*length(Controls_train)
    raw.miRNA.f.train = train_miRNA [x_train | y_train , ]
    
    
    test_miRNA  = data [ , folds$test]
    
    Controls_test = cov_pooled[folds$test, ] %>% filter ( GROUP == "Control" ) %>% rownames
    Patient_test  = cov_pooled[folds$test, ] %>% filter ( GROUP == "MDD" ) %>% rownames
    
    # filter miRs : a mir is retained only if present with more then 1
    # reads in more then 60% of either patients or controls 

    x_test = rowSums(test_miRNA[,Patient_test]  > 1) >= 0.6*length(Patient_test)
    y_test = rowSums(test_miRNA[,Controls_test] > 1) >= 0.6*length(Controls_test)
    raw.miRNA.f.test = test_miRNA [x_test | y_test, ]
    
    ### vst norm 
    dds_train_adj = DESeqDataSetFromMatrix( (raw.miRNA.f.train) , colData=cov_pooled[folds$train,] , design=~1)
    vst_train_adj = varianceStabilizingTransformation(dds_train_adj)

    dds_test_adj  = DESeqDataSetFromMatrix( (raw.miRNA.f.test) , colData=cov_pooled[folds$test,] , design=~1)
    vst_test_adj  = varianceStabilizingTransformation(dds_test_adj)
    
    ### Correct data
    if (All == TRUE){
      covmiRNA  = c( "BMI",  "Rin.miR", "AGE", "SEX" ,"Polynuclear_neutrophile", "Lymphocyte")
      
      corrected_miRNA_train = (lm(t(assay(vst_train_adj))~ Rin.miR + BMI + AGE + SEX + Polynuclear_neutrophile + Lymphocyte, 
                             data = cov_pooled[colnames(assay(vst_train_adj)) , covmiRNA]) )$residuals
      
      corrected_miRNA_test  = (lm(t(assay(vst_test_adj))~ Rin.miR + BMI + AGE + SEX + Polynuclear_neutrophile + Lymphocyte, 
                             data=cov_pooled[colnames(assay(vst_test_adj)) , covmiRNA]) )$residuals
    }else{
      covmiRNA  = c( "BMI",  "Rin.miR", "AGE" ,"Polynuclear_neutrophile", "Lymphocyte")
      
      corrected_miRNA_train = (lm(t(assay(vst_train_adj))~ Rin.miR + BMI + AGE + Polynuclear_neutrophile + Lymphocyte, 
                                  data = cov_pooled[colnames(assay(vst_train_adj)) , covmiRNA]) )$residuals
      
      corrected_miRNA_test  = (lm(t(assay(vst_test_adj))~ Rin.miR + BMI + AGE + Polynuclear_neutrophile + Lymphocyte, 
                                  data=cov_pooled[colnames(assay(vst_test_adj)) , covmiRNA]) )$residuals
    }
    return ( list ( corrected_miRNA_train = corrected_miRNA_train ,
                    corrected_miRNA_test  = corrected_miRNA_test) )
}
cv_miRNA_corr = correct_miRNA (raw.miRNA , cv_fold[[1]] , cov_pooled  )
cv_miRNA_f_corr = correct_miRNA (raw.miRNA , cv_fold_female[[1]] , cov_pooled , All = F )
cv_miRNA_m_corr = correct_miRNA (raw.miRNA , cv_fold_male[[1]] , cov_pooled , All = F )

#For correcting all folds:
# cv_miRNA_corr = lapply (cv_fold , function (x) correct_miRNA (raw.miRNA , x , cov_pooled  )  )
# cv_miRNA_f_corr = lapply (cv_fold_female , function (x) correct_miRNA (raw.miRNA , x , cov_pooled , All = FALSE )  )
# cv_miRNA_m_corr = lapply (cv_fold_male , function (x) correct_miRNA (raw.miRNA , x , cov_pooled , All = FALSE )  )

saveRDS(cv_miRNA_corr , file = "results/2_PreProcessing/cv_miRNA_corr.RDS")
saveRDS(cv_miRNA_f_corr , file = "results/2_PreProcessing/cv_miRNA_f_corr.RDS")
saveRDS(cv_miRNA_m_corr , file = "results/2_PreProcessing/cv_miRNA_m_corr.RDS")

################
###   DNAm   ###
################
myNorm.mdd = readRDS("data/myNorm.mdd.RDS") # normalised beta-values of probes
pd_mdd = readRDS("data/pd_mdd.RDS") # pd file containes metadata of samples
LeucocyteFraction.mdd = readRDS("data/LeucocyteFraction.mdd.RDS") # leucocyte fractions estimation using Houseman method

var_filter = function (DNAm_all , freq = 0.1 )  {
  ### filter according to variance 
  rvars = rowVars (DNAm_all)
  names(rvars) = rownames(DNAm_all)
  top_var = rvars %>% .[order(-.)] %>% head (n=round ( freq * length(.) ))  
  DNAm_filtered =  DNAm_all [ names(top_var) , ]
  DNAm_filtered = t(DNAm_filtered)
  return(DNAm_filtered)
}


cons_MDD_cor = function ( myNorm.mdd , pd_mdd , LeucocyteFraction.mdd, train = TRUE) {
  
  # Step one: data correction
  # using ComBat implemented in the sva package
  # the model indicating what should be "protected" variable of interest.
  #Here we do not want to protect any variable so we use a null model
  mod0 <- model.matrix(~1, data = pd_mdd)
  
  # consecutive correction removing  slide, array, sex, age, and BMI
  warning("Correcting for : Slide")
  bat <- ComBat(dat=myNorm.mdd, batch = pd_mdd$Slide, mod=mod0)
  warning("Correcting for : Array")
  if (length(table(pd_mdd$Array)) > 1){
    bat <- ComBat(dat=bat, batch=pd_mdd$Array, mod=mod0)
  }
  warning("Correcting for : Sex")
  if (length(table(pd_mdd$Sex)) > 1){
    bat <- ComBat(dat=bat, batch=pd_mdd$Sex, mod=mod0)
  }
  warning("Correcting for : Age_bin")
  bat <- ComBat(dat=bat, batch = pd_mdd$Age_bin, mod = mod0)
  warning("Correcting for : BMI.bin")
  bat <- ComBat(dat=bat, batch = pd_mdd$BMI.bin, mod = mod0)
  
  #Cell count correction
  if (train == TRUE) {
    blood = LeucocyteFraction.mdd[colnames(bat),]
    
    beta.lm = lm( t(bat) ~ CD4+CD8+MO+B+NK+GR ,
                  data=blood[bat %>% colnames , 1:6 ] )
    
    b.value.mdd.bin<- beta.lm$residuals + matrix(apply(bat, 1, mean) ,
                                                 nrow=nrow( beta.lm$residuals ), 
                                                 ncol=ncol( beta.lm$residuals ))
    b.value.mdd.bin[b.value.mdd.bin >= 1] <- 0.99999999
    b.value.mdd.bin[b.value.mdd.bin <= 0] <- 0.00000001
    correct.b.value.mdd.bin.nopreservation<-b.value.mdd.bin
    
    rownames(correct.b.value.mdd.bin.nopreservation) = pd_mdd [match(rownames(correct.b.value.mdd.bin.nopreservation) , pd_mdd$Sample_Name), "Name"]
    
    return( list ( correct.b =  correct.b.value.mdd.bin.nopreservation ,
                   model = beta.lm ) )
  }else{
    return (bat)
  }
}


correct_DNAm = function (data , pd_mdd, folds , cov_pooled, LeucocyteFraction.mdd) {
  
  warning ("Correcting train data")
  train_DNAm   = data [ , pd_mdd [ folds$train , "Sample_Name" ] ]
  corrected_DNAm_train = cons_MDD_cor (train_DNAm ,  pd_mdd [ folds$train ,  ] , LeucocyteFraction.mdd, train = TRUE)
  corrected_DNAm_train_f = var_filter ( t (corrected_DNAm_train$correct.b))
  warning ("Correcting test data")
  test_DNAm   = data [ colnames(corrected_DNAm_train_f), pd_mdd [ folds$test , "Sample_Name" ] ]
  corrected_DNAm_test = cons_MDD_cor (test_DNAm ,  pd_mdd [ folds$test ,  ] , LeucocyteFraction.mdd, train = FALSE )
  corrected_DNAm_test= corrected_DNAm_test [colnames(corrected_DNAm_train_f) , ]
  warning (paste0("Predict results for test based on coef from train: " , Sys.time() ) )
  corrected_DNAm_p =  predict (corrected_DNAm_train$model ,  cbind( corrected_DNAm_test %>% t %>% as.data.frame , LeucocyteFraction.mdd[colnames(corrected_DNAm_test), 1:6] )   )
  corrected_DNAm_p = corrected_DNAm_p [ , rownames(corrected_DNAm_test)  ]
  warning (paste0("Predict Done: " , Sys.time() ) )
  
  corrected_DNAm_test_t =  t(corrected_DNAm_test)
  residu = (corrected_DNAm_test_t- corrected_DNAm_p ) %>% as.data.frame
  b.value.mdd.bin =  as.matrix(residu)  + matrix(apply(corrected_DNAm_test_t, 2, mean) ,nrow=nrow( residu), ncol=ncol( residu ))
  
  b.value.mdd.bin[b.value.mdd.bin >= 1] <- 0.99999999
  b.value.mdd.bin[b.value.mdd.bin <= 0] <- 0.00000001
  corrected_DNAm_test_v = b.value.mdd.bin
  
  rownames(corrected_DNAm_test_v) = pd_mdd [match(rownames(corrected_DNAm_test_v) , pd_mdd$Sample_Name), "Name"] 
  corrected_DNAm_test_f = corrected_DNAm_test_v
  return ( list ( corrected_DNAm_train = corrected_DNAm_train_f ,
                  corrected_DNAm_test  = corrected_DNAm_test_f) )
}

cv_DNAm_corr = correct_DNAm (myNorm.mdd , pd_mdd, cv_fold [[1]] , cov_pooled, LeucocyteFraction.mdd)
saveRDS(cv_DNAm_corr , file = "results/2_PreProcessing/cv_DNAm_corr.RDS")

cv_DNAm_f_corr = correct_DNAm (myNorm.mdd , pd_mdd, cv_fold_female [[1]] , cov_pooled, LeucocyteFraction.mdd)
saveRDS(cv_DNAm_f_corr , file = "results/2_PreProcessing/cv_DNAm_f_corr.RDS")

cv_DNAm_m_corr = correct_DNAm (myNorm.mdd , pd_mdd, cv_fold_male [[1]] , cov_pooled, LeucocyteFraction.mdd)
saveRDS(cv_DNAm_m_corr , file = "results/2_PreProcessing/cv_DNAm_m_corr.RDS")
