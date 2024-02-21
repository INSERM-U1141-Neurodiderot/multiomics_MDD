library(dplyr)
library(stringr)
library(DESeq2)
library(ChAMP)
library(sva)
library(reticulate)

### Set your working directory
setwd("/env/export/v_cng_n07_scratch/v_scratch_math_stats/gloaguen/neurodiderot/repro_paper/multiomics_MDD/3_multiomics/")

### input: covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

### input: cv folds
cv_fold = readRDS("results/1_CrossValidation/cv_fold.RDS")

################
###   mRNA   ###
################

data = readRDS('data/data_mRNA.RDS')

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

#For correcting all folds:
cv_mRNA_corr = lapply (cv_fold , function (x) correct_mRNA (data , x , cov_pooled  )  )

saveRDS (cv_mRNA_corr , file = "results/2_PreProcessing/cv_mRNA_corr.RDS")

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
    col = colnames(corrected_miRNA_test)[colnames(corrected_miRNA_test) %in% colnames(corrected_miRNA_train)]
    corrected_miRNA_test = corrected_miRNA_test[, col]
    corrected_miRNA_train = corrected_miRNA_train[, col]
    return ( list ( corrected_miRNA_train = corrected_miRNA_train ,
                    corrected_miRNA_test  = corrected_miRNA_test) )
}

#For correcting all folds:
cv_miRNA_corr = lapply (cv_fold , function (x) correct_miRNA (raw.miRNA , x , cov_pooled  )  )

saveRDS(cv_miRNA_corr , file = "results/2_PreProcessing/cv_miRNA_corr.RDS")

################
###   DNAm   ###
################
var_filter = function (DNAm_all , freq = 0.1 )  {
  ### filter according to variance 
  rvars = rowVars (DNAm_all)
  names(rvars) = rownames(DNAm_all)
  top_var = rvars %>% .[order(-.)] %>% head (n=round ( freq * length(.) ))  
  DNAm_filtered =  DNAm_all [ names(top_var) , ]
  DNAm_filtered = t(DNAm_filtered)
  return(DNAm_filtered)
}

reticulate::virtualenv_create(envname = "repro_multiomics")
np = import("numpy")
CombatModel = import_from_path("neurocombat_sklearn", path = "/env/export/v_home/q_unix/agloague/.virtualenvs/repro_multiomics/lib/python3.8/site-packages")
neurocombat_transfert = function(model , DataTrain , DataTest , CovData , Cov , Train , Test) {
  Model  = model 
  
  if (length(Cov) > 1 ) {    
    
    X_train = Model$fit_transform( DataTrain    ,
                                   array_reshape(CovData [ Train ,  Cov [1] ] , c(-1, 1)) ,
                                   array_reshape(CovData [ Train ,   Cov [2] ] , c(-1, 1)) )
    
    X_test = Model$transform(      DataTest   ,
                                   array_reshape(CovData [ Test ,  Cov [1] ] , c(-1, 1)) ,
                                   array_reshape(CovData [ Test ,   Cov [2] ] , c(-1, 1)) )
  }
  if (length(Cov) == 1 ) {    
    
    X_train = Model$fit_transform(DataTrain    ,
                                  as.data.frame(as.factor(CovData [ Train ,  Cov ])))
    
    X_test = Model$transform(DataTest   ,
                             as.data.frame(as.factor(CovData [ Test ,  Cov ]))) 
  }
  
  return( list (X_train = X_train , 
                X_test = X_test ) )
  
}

neurocombat_correct = function(model , Data  , CovData , Cov ) {
  Model  = model 
  
  Data_corrected = Model$fit_transform(Data ,
                                       as.data.frame(as.factor(CovData[ , Cov])))
  
  return( Data_corrected )
  
}

correction_DNAm = function(cv_fold, DNAm.npy, i, pd_mdd2, LeucocyteFraction.mdd, freq, CovDNAm, CovSep){
        DNAm.npy = var_filter(DNAm.npy, freq)
        
        # Creating models
        model  = CombatModel$CombatModel()
        modelFtrain = CombatModel$CombatModel()
        modelFetest = CombatModel$CombatModel()

        Combat_DNAm = list (X_train =  DNAm.npy [cv_fold[[i]]$train , ]  ,  X_test  =  DNAm.npy [cv_fold[[i]]$test ,  ] ) 
        for (Cov in CovDNAm) {
                model  = CombatModel$CombatModel()
        
                Combat_DNAm = neurocombat_transfert ( model =  model ,  DataTrain =  Combat_DNAm$X_train  , 
                                          DataTest  =  Combat_DNAm$X_test ,  CovData = pd_mdd2  ,  Cov =  Cov , 
                                          Train = cv_fold[[i]]$train , 
                                          Test  = cv_fold[[i]]$test  )
        }
        
        modelFtrain = CombatModel$CombatModel()
        modelFetest = CombatModel$CombatModel()
        
        for (Cov in CovSep){
                DNAm_Train_cor = neurocombat_correct (model = modelFtrain, Data = Combat_DNAm$X_train , CovData = pd_mdd2 [cv_fold[[i]]$train ,], Cov = Cov)
                DNAm_Test_cor = neurocombat_correct (model = modelFetest, Data = Combat_DNAm$X_test , CovData = pd_mdd2 [cv_fold[[i]]$test ,], Cov = Cov)
        }
        DNAm_Train = as.data.frame (DNAm_Train_cor)
        DNAm_Test = as.data.frame (DNAm_Test_cor)
        
        colnames(DNAm_Train) = colnames(DNAm_Test) = colnames(DNAm.npy)
        rownames (DNAm_Train) = cv_fold[[i]]$train
        rownames (DNAm_Test) = cv_fold[[i]]$test
        
        b.value.mdd.bin = NULL
        beta.lm = lm (  as.matrix(DNAm_Train)  ~ CD4+CD8+MO+B+NK+GR , data = LeucocyteFraction.mdd[ rownames(DNAm_Train)  , 1:6 ]  )
        b.value.mdd.bin<- beta.lm$residuals + matrix(apply(DNAm_Train, 2, mean) ,nrow=nrow( beta.lm$residuals ), ncol=ncol( beta.lm$residuals ))
        
        b.value.mdd.bin[b.value.mdd.bin >= 1] <- 0.99999999
        b.value.mdd.bin[b.value.mdd.bin <= 0] <- 0.00000001
        
        DNAm_Train_c = b.value.mdd.bin
        
        b.value.mdd.res= NULL
        corrected_DNAm_p =  predict ( beta.lm ,  cbind( DNAm_Test %>% as.data.frame , LeucocyteFraction.mdd[rownames(DNAm_Test), 1:6]  )   )
        residu = (as.matrix(DNAm_Test) - corrected_DNAm_p) %>% as.data.frame
        b.value.mdd.res<- as.matrix(residu) + matrix(apply(DNAm_Test, 2, mean) ,nrow=nrow( residu), ncol=ncol( residu )) 
        b.value.mdd.res[b.value.mdd.res >= 1] <- 0.99999999
        b.value.mdd.res[b.value.mdd.res <= 0] <- 0.00000001
        
        return ( list ( corrected_DNAm_train = DNAm_Train_c ,
                        corrected_DNAm_test  = b.value.mdd.res) )
}

myNorm.mdd            = readRDS("../../../multiomics_MDD-main/3_multiomics/data/myNorm.mdd.RDS") # normalised beta-values of probes
pd_mdd                = readRDS("data/pd_mdd.RDS") # pd file containes metadata of samples
LeucocyteFraction.mdd = readRDS("data/LeucocyteFraction.mdd.RDS") # leucocyte fractions estimation using Houseman method

colnames(myNorm.mdd) = rownames(pd_mdd)[match(colnames(myNorm.mdd), pd_mdd$Sample_Name)]
myNorm.mdd           = myNorm.mdd[, -which(is.na(colnames(myNorm.mdd)))]

tmp_rownames                    = rownames(pd_mdd)[match(rownames(LeucocyteFraction.mdd), pd_mdd$Sample_Name)]
idx_to_rm                       = which(is.na(tmp_rownames))
LeucocyteFraction.mdd           = LeucocyteFraction.mdd[-idx_to_rm, ]
rownames(LeucocyteFraction.mdd) = tmp_rownames[-idx_to_rm]

CovDNAm      = c("Array" ,  "Sex" , "BMI.bin" , "Age_bin"  )
CovSep       = c("Slide")
cv_DNAm_corr = lapply (1:length(cv_fold), function(x) correction_DNAm(cv_fold      = cv_fold, 
                                                             DNAm.npy              = myNorm.mdd, 
                                                             i                     = x, 
                                                             pd_mdd2               = pd_mdd, 
                                                             LeucocyteFraction.mdd = LeucocyteFraction.mdd, 
                                                             freq                  = 0.1,
                                                             CovDNAm               = CovDNAm,
							     CovSep                = CovSep))

saveRDS(cv_DNAm_corr , file = "results/2_PreProcessing/cv_DNAm_corr_filter_10.RDS")
