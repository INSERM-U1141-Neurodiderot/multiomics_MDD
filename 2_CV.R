#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
#########################
np = import("numpy")
CombatModel = import("neurocombat_sklearn")
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
    
    X_train = Model$fit_transform( DataTrain    ,
                                   array_reshape(CovData [ Train ,  Cov ] , c(-1, 1)) 
    )
    
    X_test = Model$transform(      DataTest   ,
                                   array_reshape(CovData [ Test ,  Cov  ] , c(-1, 1)) 
    )
  }
  
  return( list (X_train = X_train , 
                X_test = X_test ) )
  
}

neurocombat_correct = function(model , Data  , CovData , Cov ) {
  Model  = model 
  
  Data_corrected = Model$fit_transform( Data ,
                                        array_reshape( CovData [  ,   Cov  ] , c(-1, 1)) )
  
  return( Data_corrected )
  
}

CovDNAm = c("Array" ,  "Sex" , "BMI.bin" , "Age_bin"  )
blood = readRDS("data/LeucocyteFraction.mdd_2.RDS")
diff_DNAm = function (DNAm_T , cov_pooled) {
  
  # Step four: DMP analysis
  # This will use limma package: there for it is better to transform beta-values to M-values (Du et al. 2010)    
  corrected_m.value.mdd.bin.nopreservation = lumi::beta2m(DNAm_T )
  
  #DMP analysis
  # input data: pd file, m-values
  # compare.group: we compare mdd group vs control group
  myDMP.mdd.bin.nopreservation = champ.DMP( beta         = t(corrected_m.value.mdd.bin.nopreservation) , 
                                            pheno         =  cov_pooled [rownames(corrected_m.value.mdd.bin.nopreservation) , "GROUP" ] %>% as.character , 
                                            #  compare.group = c("Control","MDD"), 
                                            adjPVal       = 1, 
                                            adjust.method = "BH", 
                                            arraytype     = "EPIC")
  
  return(myDMP.mdd.bin.nopreservation [[1]] )
  
}

for (i in 1:length(cv_fold) ) {
  
  # Creating models
  model  = CombatModel$CombatModel()
  modelFtrain = CombatModel$CombatModel()
  modelFetest = CombatModel$CombatModel()
  
  Combat_DNAm = list (X_train =  DNAm.npy [cv_fold[[i]]$train , ]  , 
                      X_test  =  DNAm.npy [cv_fold[[i]]$test ,  ] ) 
  for (Cov in CovDNAm) {
    model  = CombatModel$CombatModel()
    
    Combat_DNAm = neurocombat_transfert ( model =  model , 
                                          DataTrain =  Combat_DNAm$X_train  , 
                                          DataTest  =  Combat_DNAm$X_test , 
                                          CovData = pd_mdd2  , 
                                          Cov =  Cov , 
                                          Train = cv_fold[[i]]$train , 
                                          Test  = cv_fold[[i]]$test  )
  }
  
  modelFtrain = CombatModel$CombatModel()
  modelFetest = CombatModel$CombatModel()
  
  DNAm_Train_cor = neurocombat_correct (model = modelFtrain,
                                        Data = Combat_DNAm$X_train ,
                                        CovData = pd_mdd2 [cv_fold[[i]]$train ,],
                                        Cov = "Slide")
  
  
  DNAm_Test_cor = neurocombat_correct (model = modelFetest,
                                       Data = Combat_DNAm$X_test ,
                                       CovData = pd_mdd2 [cv_fold[[i]]$test ,],
                                       Cov = "Slide")
  
  
  
  
  
  DNAm_Train = as.data.frame (DNAm_Train_cor)
  DNAm_Test = as.data.frame (DNAm_Test_cor)
  
  
  colnames(DNAm_Train) = colnames(DNAm_Test) = colnames(DNAm.npy)
  rownames (DNAm_Train) = cv_fold[[i]]$train
  rownames (DNAm_Test) = cv_fold[[i]]$test
  
  
  b.value.mdd.bin = NULL
  beta.lm = lm (  as.matrix(DNAm_Train)  ~ CD4+CD8+MO+B+NK+GR , data = blood[ rownames(DNAm_Train)  , 1:6 ]  )
  b.value.mdd.bin<- beta.lm$residuals + matrix(apply(DNAm_Train, 2, mean) ,nrow=nrow( beta.lm$residuals ), ncol=ncol( beta.lm$residuals ))
  
  b.value.mdd.bin[b.value.mdd.bin >= 1] <- 0.99999999
  b.value.mdd.bin[b.value.mdd.bin <= 0] <- 0.00000001
  
  DNAm_Train_c = b.value.mdd.bin
  
  b.value.mdd.res= NULL
  corrected_DNAm_p =  predict ( beta.lm ,  cbind( DNAm_Test %>% as.data.frame , blood[rownames(DNAm_Test), 1:6]  )   )
  residu = (as.matrix(DNAm_Test) - corrected_DNAm_p) %>% as.data.frame
  b.value.mdd.res<- as.matrix(residu) + matrix(apply(DNAm_Test, 2, mean) ,nrow=nrow( residu), ncol=ncol( residu )) 
  b.value.mdd.res[b.value.mdd.res >= 1] <- 0.99999999
  b.value.mdd.res[b.value.mdd.res <= 0] <- 0.00000001
  
  DNAm_Test_c = b.value.mdd.res
  
  
  Diff_Train =  diff_DNAm  (DNAm_Train_c , cov_pooled)
  
  saveRDS (DNAm_Train_c , file = paste0("results/CV/Diff/" , train_nms  [[i]]) )
  saveRDS (DNAm_Test_c , file = paste0("results/CV/Diff/" , test_nms [[i]]) )
  saveRDS (Diff_Train ,  file = paste0("results/CV/Diff/" , DNAmDiffCv_nms [[i]]))
}
#################################################################

DNAmDif = lapply (DNAmDiffCv_nms , function(x) readRDS (paste0("results/CV/Diff/" , x )) )

allCG = lapply (1: length (DNAmDif)  , function (x) DNAmDif [[x]] %>% filter (P.Value < 10^-4) %>% row.names) %>% unlist 

TrainDNAm = list ()
TestDNAm = list ()

for (i in 1:length(DNAmDiffCv_nms)) {
  diff = readRDS (paste0("results/CV/Diff/" , DNAmDiffCv_nms [[i]]  ) )
  DNAm_train =  readRDS( paste0("results/CV/Diff/" , train_nms  [[i]]) )
  DNAm_test =  readRDS(  paste0("results/CV/Diff/" , test_nms [[i]]) ) 
  
  TrainDNAm = append (TrainDNAm , list(DNAm_train [,diff%>% filter (P.Value < 10^-4) %>% row.names ] ) )
  TestDNAm =  append (TestDNAm , list(DNAm_test [, diff %>% filter (P.Value < 10^-4) %>% row.names ] ) )
}

saveRDS(TrainDNAm , "results/CV/Diff/TrainDNAm0.0001.RDS" )
saveRDS(TestDNAm  , "results/CV/Diff/TestDNAm0.0001.RDS" )

saveRDS(TrainDNAm , "results/CV/Diff/TrainDNAm0.01.RDS" )
saveRDS(TestDNAm  , "results/CV/Diff/TestDNAm0.01.RDS" )

##########################################################################################################################

cv_miRNA_corr = readRDS("results/CV/cv_miRNA_corr.RDS")
cv_mRNA_corr = readRDS("results/CV/cv_mRNA_corr.RDS")

trainDiff_cvs_dataset = lapply (1:25 , function(x) list ( list (  list( DNAm = TrainDNAm      [[x]]   ) )  ,
                                                          list (  list( mRNA = cv_mRNA_corr   [[x]][[1]] ) ),
                                                          list ( list(  mRNA = cv_mRNA_corr   [[x]][[1]] ), list( DNAm = TrainDNAm [[x]] )),
                                                          list ( list(  miRNA = cv_miRNA_corr [[x]][[1]] ) ) ,
                                                          list (  list( miRNA = cv_miRNA_corr [[x]][[1]]),  list( DNAm = TrainDNAm [[x]])),
                                                          list (  list( miRNA = cv_miRNA_corr [[x]][[1]] ),list(  mRNA = cv_mRNA_corr [[x]][[1]] )),
                                                          list ( list(  miRNA = cv_miRNA_corr [[x]][[1]]) , list(mRNA = cv_mRNA_corr [[x]][[1]]) , list( DNAm = TrainDNAm [[x]]) ) 
) )


testDiff_cvs_dataset = lapply (1:25 , function(x) list ( list (   list(DNAm = TestDNAm      [[x]]) ) ,
                                                         list (  list( mRNA = cv_mRNA_corr   [[x]][[2]] )),
                                                         list (  list( mRNA = cv_mRNA_corr   [[x]][[2]]) , list( DNAm = TestDNAm [[x]] )),
                                                         list (  list( miRNA = cv_miRNA_corr [[x]][[2]] ) ),
                                                         list (  list( miRNA = cv_miRNA_corr [[x]][[2]]), list( DNAm = TestDNAm [[x]])),
                                                         list (  list( miRNA = cv_miRNA_corr [[x]][[2]] ), list(mRNA = cv_mRNA_corr [[x]][[2]] )),
                                                         list (  list( miRNA = cv_miRNA_corr [[x]][[2]]), list(mRNA = cv_mRNA_corr [[x]][[2]]) ,  list( DNAm = TestDNAm [[x]]) ) 
) ) 


save (trainDiff_cvs_dataset , testDiff_cvs_dataset , file = "DiffDataP0.0001.RData")

####################################
###   Create CrossVal datasets   ###
####################################

cv_DNAm_corr_train = lapply (paste0("results/CV/" , train_nms ) , readRDS) 
cv_DNAm_corr_test = lapply (paste0("results/CV/" , test_nms ) , readRDS) 

cv_DNAm =  list.files(path= "results/CV", pattern="^Cv", recursive  = TRUE,
                      all.files  =TRUE, full.names =TRUE) %>% . [grep ("*test*" , .)]

cv_DNAm_corr_test = lapply (cv_DNAm , readRDS)

cv_miRNA_corr = readRDS("results/CV/cv_miRNA_corr.RDS")
cv_mRNA_corr = readRDS("results/CV/cv_mRNA_corr.RDS")

Hamilton = cov_pooled %>% dplyr::select(starts_with('Hamilton'))

cv_all_train = lapply (1:25 , function (x)  list (miRNA = cv_miRNA_corr[[x]] [[1]] ,
                                                  mRNA  = cv_mRNA_corr[[x]] [[1]] ,
                                                  DNAm  = cv_DNAm_corr_train [[x]],
                                                  Cov   = Hamilton [cv_fold[[x]][[1]], ]) 
)


cv_all_test = lapply (1:25 , function (x)   list (miRNA = cv_miRNA_corr[[x]] [[2]] ,
                                                  mRNA  = cv_mRNA_corr[[x]] [[2]] ,
                                                  DNAm  = cv_DNAm_corr_test [[x]],
                                                  Cov   = Hamilton [cv_fold[[x]][[2]], ]) 
)

saveRDS (cv_all_train , "results/CV/cv_all_train.RDS")
saveRDS (cv_all_test , "results/CV/cv_all_test.RDS")

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 1:5 , function (x){
  saveRDS (cv_all_train [[x]] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(3,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,2,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,3,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,3,4)] , paste0( "results/CV/Allcv1P" , x , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 6:10 , function (x){
  saveRDS (cv_all_train [[x]] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(3,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,2,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,3,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,3,4)] , paste0( "results/CV/Allcv2P" , (x -5) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 11:15 , function (x){
  saveRDS (cv_all_train [[x]] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(3,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,2,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,3,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,3,4)] , paste0( "results/CV/Allcv3P" , (x -10) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 16:20 , function (x){
  saveRDS (cv_all_train [[x]] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(3,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,2,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,3,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,3,4)] , paste0( "results/CV/Allcv4P" , (x -15) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 21:25 , function (x){
  saveRDS (cv_all_train [[x]] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(3,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,2,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_train [[x]] [c(1,3,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_train [[x]] [c(2,3,4)] , paste0( "results/CV/Allcv5P" , (x -20) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 1:5 , function (x){
  saveRDS (cv_all_test [[x]] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(3,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,2,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,3,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,3,4)] , paste0( "results/TEST/Allcv1P" , x , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 6:10 , function (x){
  saveRDS (cv_all_test [[x]] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(3,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,2,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,3,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,3,4)] , paste0( "results/TEST/Allcv2P" , (x -5) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 11:15 , function (x){
  saveRDS (cv_all_test [[x]] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(3,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,2,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,3,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,3,4)] , paste0( "results/TEST/Allcv3P" , (x -10) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 16:20 , function (x){
  saveRDS (cv_all_test [[x]] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(3,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,2,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,3,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,3,4)] , paste0( "results/TEST/Allcv4P" , (x -15) , "_MDD_omic_mRNADNAm.RDS")) } )

#Female_MDD_omic_miRNAmRNADNAm
lapply ( 21:25 , function (x){
  saveRDS (cv_all_test [[x]] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_miRNAmRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_miRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_mRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(3,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_DNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,2,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_miRNAmRNA.RDS"))
  saveRDS (cv_all_test [[x]] [c(1,3,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_miRNADNAm.RDS"))
  saveRDS (cv_all_test [[x]] [c(2,3,4)] , paste0( "results/TEST/Allcv5P" , (x -20) , "_MDD_omic_mRNADNAm.RDS")) } )

#################################################################
nemo_f = function (data, tl2 = tl , names = c("miRNA", "mRNA", "DNAm", "AllOmics")) 
{
  setClass(Class = "NemoRes", representation(graph = "matrix", 
                                             clustering = "integer"))
  nemo.clust = function(omics.list, num.clusters = NULL, num.neighbors = NA) {
    if (is.null(num.clusters)) {
      num.clusters = NA
    }
    graph = nemo.affinity.graph(omics.list, k = num.neighbors)
    if (is.na(num.clusters)) {
      num.clusters = nemo.num.clusters(graph)
    }
    clustering = spectralClustering(graph, num.clusters)
    names(clustering) = colnames(graph)
    return(new("NemoRes", graph = graph, clustering = clustering))
  }
  nemo.clust_momix = lapply(data, function(x) nemo.clust(list(x), 
                                                         num.clusters = 2, num.neighbors = NA))
  nemo_all = nemo.clust(data, num.clusters = 2, num.neighbors = NA)
  nemo_group = do.call(cbind, lapply(nemo.clust_momix, function(x) {
    x@clustering
  })) %>% as.data.frame
  nemo_group$factALL = nemo_all@clustering
  nemo_group$TrueLabel = tl2
  ARIs_nemo = lapply(nemo_group[, 1:(ncol(nemo_group) - 1)], 
                     function(x) round (adjustedRandIndex(x, nemo_group[, ncol(nemo_group)]) ,4)   )
  names(ARIs_nemo) = names
  return(list(ARIs = ARIs_nemo, Clusters = nemo_group))
}

SNF_f = function (Data,
                  K = 30 ,
                  alpha = 0.6 ,
                  T = 50  ,
                  C = 2  ,
                  tl2 = tl , 
                  color = c("#ffe935","#e289bf"),
                  names = c("miRNA", "mRNA", "DNAm", "AllOmics")) {
  Dist.MOMIX_metaG = lapply(lapply(Data, t), function(x) {
    (dist2((as.matrix(x)), as.matrix(x)))^(1/2)
  })
  W.MOMIX_metaG = lapply(Dist.MOMIX_metaG, function(x) {
    affinityMatrix(x, K, alpha)
  })
  W_all_MOMIX_metaG = SNF(W.MOMIX_metaG, K, T)
  group_mom = lapply(W.MOMIX_metaG, function(x) {
    spectralClustering(x, C)
  })
  if (length(Data) > 1) {
    group_mom = append(group_mom, list(factALL = spectralClustering(W_all_MOMIX_metaG, 
                                                                    C)))
  }
  levels(tl2) = c(1, 2)
  M_label_ = t(do.call(rbind, group_mom)) %>% as.data.frame
  M_label_$TrueLabel = tl2
  rownames(M_label_) = rownames(W.MOMIX_metaG[[1]])
  if ( length(Data) == 1) {
    M_label_$factALL = M_label_ [ , 1 ]
  }
  if (length(Data) > 1) {
    ARIs = lapply(M_label_[, 1:(ncol(M_label_) - 1)], function(x) adjustedRandIndex(x, 
                                                                                    M_label_[, ncol(M_label_)]))
    names(ARIs) = names          
    
  }
  else {
    ARIs = adjustedRandIndex(M_label_[, 1], M_label_[, 2])
    names(ARIs) = paste0(names[1:length(Data)] , "_comb")
  }
  
  return(list(ARIs = ARIs, Clusters = M_label_))
}
SNF_f = function (Data,
                  K = 30 ,
                  alpha = 0.6 ,
                  T = 50  ,
                  C = 2  ,
                  tl2 = tl , 
                  color = c("#ffe935","#e289bf"),
                  names = c("miRNA", "mRNA", "DNAm", "AllOmics")) {
  Dist.MOMIX_metaG = lapply(lapply(Data, t), function(x) {
    (dist2((as.matrix(x)), as.matrix(x)))^(1/2)
  })
  W.MOMIX_metaG = lapply(Dist.MOMIX_metaG, function(x) {
    affinityMatrix(x, K, alpha)
  })
  W_all_MOMIX_metaG = SNF(W.MOMIX_metaG, K, T)
  group_mom = lapply(W.MOMIX_metaG, function(x) {
    spectralClustering(x, C)
  })
  if (length(Data) > 1) {
    group_mom = append(group_mom, list(factALL = spectralClustering(W_all_MOMIX_metaG, 
                                                                    C)))
  }
  levels(tl2) = c(1, 2)
  M_label_ = t(do.call(rbind, group_mom)) %>% as.data.frame
  M_label_$TrueLabel = tl2
  rownames(M_label_) = rownames(W.MOMIX_metaG[[1]])
  if ( length(Data) == 1) {
    M_label_$factALL = M_label_ [ , 1 ]
  }
  if (length(Data) > 1) {
    ARIs = lapply(M_label_[, 1:(ncol(M_label_) - 1)], function(x) adjustedRandIndex(x, 
                                                                                    M_label_[, ncol(M_label_)]))
    names(ARIs) = names          
    
  }
  else {
    ARIs = adjustedRandIndex(M_label_[, 1], M_label_[, 2])
    names(ARIs) = paste0(names[1:length(Data)] , "_comb")
  }
  
  return(list(ARIs = ARIs, Clusters = M_label_))
}
generate_omics_dataset = function (omics , dat_list , x) {
  
  omics_data_mcia   = omics_data_sets (omics , dat_list [[x]] $top_mcia )
  omics_data_mofa   = omics_data_sets (omics , dat_list [[x]] $top_mofa )
  omics_data_intNMF = omics_data_sets (omics , dat_list [[x]] $top_intNMF )
  omics_data_jive   = omics_data_sets (omics , dat_list [[x]] $top_jive )
  omics_data_scikit = omics_data_sets (omics , dat_list [[x]] $top_scikit )
  omics_data_rgcca  = omics_data_sets (omics , dat_list [[x]] $top_rgcca )
  
  
  omics_datasets = list (mcia = omics_data_mcia,
                         mofa = omics_data_mofa,
                         intNMF = omics_data_intNMF ,
                         jive = omics_data_jive,
                         scikit = omics_data_scikit,
                         rgcca = omics_data_rgcca)
  
  return(omics_datasets)
}
df_top_features =  list.files(path= "01_Results/CV_train_results", 
                              pattern=".RDS",
                              recursive = TRUE,
                              all.files=TRUE,
                              full.names=TRUE)


### Read the data
dat_list = list()

for (i in 1:length(df_top_features)) {
  F = readRDS (df_top_features[i])
  dat_list [[i]]  =   F
}


df_to_names =  list.files(path= "01_Results/CV_train_results", 
                          pattern=".RDS",
                          recursive = TRUE,
                          all.files=TRUE,
                          full.names=FALSE)


df_to_names= lapply (df_to_names , function (x) gsub ("_.+" , "",  x) )


names(dat_list) = df_to_names


Allcv_train = lapply (df_to_names , function (x) readRDS(  paste0( "00_Data/CV/" , x , "_MDD_omic_miRNAmRNADNAm.RDS")) )
Allcv_test  =lapply (df_to_names , function (x)  readRDS(  paste0( "00_Data/TEST/" , x , "_MDD_omic_miRNAmRNADNAm.RDS")) )





cg_cv = lapply ( Allcv_train , function (x) colnames(x[[3]]) )


cg_cv_unique = cg_cv %>% unlist %>% unique


cg_cv_c = lapply (cg_cv , function (x) cg_cv_unique %in% x %>% as.double  )


cg_cv_cdf = cg_cv_c %>% do.call(cbind , .)
rownames(cg_cv_cdf) = cg_cv_unique


cg_cv_p = (cg_cv_cdf %>% rowSums * 100) / 25


cg_cv_p %>% as.data.frame %>%
  ggplot (. , aes(x = .)) + geom_histogram()+
  ggtitle ("DNAm before feature selection across the 25 split of CV ")



names(Allcv_test ) = df_to_names
names(Allcv_train ) = df_to_names




# construct omics datastest from the top features 
omics_data_sets = function (omics , top_meth ) {
  k = list()
  omics_data   = mapply ( function (x , y)  { lapply (x , function(z)  { k [[z]] = omics [[z]] [,  intersect( colnames(omics [[z]] ) , top_meth [[y]] [[z]] )  ]   ; k } )   }   , x = lapply ( top_meth , names )  , y = 1:length(top_meth) )
  return(omics_data)
}




train_cvs_dataset [[1]]$Diff = trainDiff_cvs_dataset[[1]]


train_cvs_dataset = lapply (1:25 , function (x) generate_omics_dataset (Allcv_train [[x]], dat_list , x ) )
test_cvs_dataset  = lapply (1:25 , function (x) generate_omics_dataset (Allcv_test [[x]] , dat_list , x ) )


lists_names = lapply (1:25 , function (x) lapply (dat_list [[x]] $top_mcia , names) )
Group_a = lapply (1:25 , function (x) cov_pooled [Allcv_test[[x]]$miRNA %>% rownames , "GROUP" ] %>% as.character )


input_names = lapply (1:25 , function (x) paste0("input_" , lapply(lists_names [[x]] , function (x) paste0(x , collapse = "_")) ) ) 



names(Group_a) = names(Allcv_test)


load ("DiffDataP0.0001.RData")


### Add Diff dataset 
for (i in 1:25) {train_cvs_dataset [[i]]$Diff =  trainDiff_cvs_dataset[[i]]  }
for (i in 1:25) {test_cvs_dataset [[i]]$Diff  = testDiff_cvs_dataset[[i]]    }



#### Aris results from clusturing : NEMO
input_names = lapply (1:25 , function (x) paste0("input_" , lapply(lists_names [[x]] , function (x) paste0(x , collapse = "_")) ) ) 
nemo_res    = lapply (1:25 ,  function (G) lapply (test_cvs_dataset [[G]] , function (Z) lapply ( Z , function (Y) { nemo_f (   lapply (Y  ,function(x)  t(x [[1]] )  )  , tl2 = Group_a [[G]] , names =  c( unlist ( lapply (Y ,function(x)  names (x  )  ) ) , paste0( c(unlist ( lapply (Y ,function(x)  names (x  )  ) ) , "comb"), collapse = "_" ) ) )}   ) ) )
nemo_res    = lapply (1:25 ,  function (Y) lapply (nemo_res [[Y]] , function (x) {names(x) =  input_names [[Y]] ; x} ) )

#### Aris results from clusturing : SNF
snf_res = lapply (1:25 ,  function (G) lapply (test_cvs_dataset [[G]], function (Z) lapply ( Z , function (Y) { SNF_f (   lapply (Y  ,function(x)  t(x [[1]] )  )  , tl2 = Group_a [[G]] , names =  c( unlist ( lapply (Y ,function(x)  names (x  )  ) ) , paste0( c(unlist ( lapply (Y ,function(x)  names (x  )  ) ) , "comb"), collapse = "_" )  )   )}   ) ) )
snf_res = lapply (1:25 ,  function (Y) lapply (snf_res [[Y]] , function (x) {names(x) =  input_names [[Y]] ; x} ) )

