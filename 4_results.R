SNF_Fun = function (Data ,params) {
  # Normalization 
  NomalizedOm = lapply (Data , standardNormalization ) 
  # Compute Distance matrix
  DistM  = lapply( NomalizedOm  , function(x) { ( dist2 (  as.matrix(x) , as.matrix(x) ) ) ^(1/2) } )
  # Compute affinity Matrix
  W.M = lapply( DistM, function(x) { affinityMatrix (x, paramsalpha) } )
  # Apply SNF
  W  =  SNF(W.M, K = paramst)
  #  Spectral Clustering
  G = spectralClustering( W , K = params$Nclust)
  return(G)
}

ConfM = function (predLab , TrueLabs , name) { 
  # caret::confusionMatrix
  ConfMTest =  confusionMatrix( as.factor (TrueLabs) , as.factor (predLab)  )
  byClass   =  ConfMTestoverall
  # rename Results for display 
  names (byClass) =  paste0(  paste0( name , "class") ,   lapply (names (byClass) ,  function (x) gsub(" " , "_" , x ) ) )
  names (overall) =  paste0(  paste0( name , "overall") , lapply (names (overall) ,  function (x) gsub(" " , "_" , x ) ) )
  return(  c (byClass , overall ) )
}

SNF_Pred = function (DataTrain , DataTest , labTrain , labTest , params , TrueLabs = FALSE ){
  
  # if single omic => duplicate 
  if ( length( DataTrain ) == 1 ) {
    DataTrain = c( DataTrain, DataTrain )
    DataTest  = c( DataTest, DataTest )
  }
  ### Use the real Labs 
  if (TrueLabs == TRUE ) {
    TrainLab = labTrain %>% as.factor %>% as.double
  }else{
    # Label from  train 
    TrainLab = SNF_Fun ( DataTrain , params  )
  }
  
  # predict labels based on SNF laabs of train 
  newLabel =     groupPredict( train =   DataTrain , 
                               test  =   DataTest  ,
                               TrainLab  ,  
                               K = paramsalpha ,
                               t = params$t ) 
  
  # evalutae labs 
  TrueLabTrain = labTrain %>% as.factor %>% as.double 
  TrueLabTest  = labTest %>% as.factor %>% as.double 
  
  # get the predicted Labs for Test                                      
  newTestLabel = newLabel [( length ( TrueLabTrain ) + 1 ):length ( newLabel ) ] %>% as.factor
  levels(newTestLabel) = c('1' ,'2')
  
  # Confusion Matrix 
  ConfMTest  =  ConfM( newTestLabel , TrueLabTest  , "test_"  )
  ConfMTrain =  ConfM( TrainLab     , TrueLabTrain , "train_"  )        
  
  
  AriTrain = adjustedRandIndex (TrueLabTrain , TrainLab  )
  AriTest  = adjustedRandIndex (TrueLabTest  , newTestLabel  )
  names(AriTrain) = "AriTrain"
  names(AriTest)  = "AriTest"
  return ( c( ConfMTest , ConfMTrain , AriTrain , AriTest )  )
}



SNFres = function (CV,MET,OM , labTrainS, labTestS , params)  {
  
  warning ( paste0 ("CV: " , CV , "MET: " , MET , "OM: " , OM ))
  DataTrain = lapply  ( train_cvs_dataset [[CV]] [[MET]] [[OM]]  , function (x) x [[1]] )
  DataTest  = lapply  ( test_cvs_dataset [[CV]] [[MET]] [[OM]]   , function (x) x [[1]] )
  
  labTrain =  labTrainS  [[CV]] 
  labTest  =  labTestS   [[CV]]
  
  SNFPred = SNF_Pred (DataTrain , DataTest , labTrain , labTest )                  
  return ( append ( list (cv = CV  ,
                          method = names( train_cvs_dataset [[CV]] ) [[MET]]  ,
                          omics = paste0 (lapply ( test_cvs_dataset [[CV]] [[MET]] [[OM]]  , names  )  , collapse = "_" ) )  ,
                    SNFPred ) 
  )
  
}


labTrainS = lapply (1:25 , function (x) cov_pooled [Allcv_train[[x]]$miRNA %>% rownames , "GROUP" ] %>% as.character )
labTestS =  lapply (1:25 , function (x) cov_pooled [Allcv_test[[x]]$miRNA %>% rownames , "GROUP" ] %>% as.character )


params = list ( K = 60 ,
                alpha = 0.6 ,
                t = 50 ,
                Nclust  = 2)






NoFselect = list()

t = lapply ( 1:25 , function (x) {  NoFselect [x] = list ( list ( list ( list (DNAm = Allcv_train [[x]] [["DNAm"]] )   ) ,
                                                                  list ( list ( mRNA = Allcv_train [[x]] [["mRNA"]]) )   ,
                                                                  list ( list ( mRNA = Allcv_train [[x]] [["mRNA"]]) ,
                                                                         list ( DNAm = Allcv_train [[x]] [["DNAm"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_train [[x]] [["miRNA"]]) )   ,
                                                                  list ( list ( miRNA = Allcv_train [[x]] [["miRNA"]]) ,
                                                                         list ( DNAm = Allcv_train [[x]] [["DNAm"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_train [[x]] [["miRNA"]]) ,
                                                                         list ( mRNA = Allcv_train [[x]] [["mRNA"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_train [[x]] [["miRNA"]]) ,
                                                                         list ( mRNA = Allcv_train [[x]] [["mRNA"]] )   , 
                                                                         list ( DNAm = Allcv_train [[x]] [["DNAm"]]) )  )  )   } )


NoFselect = list()

s = lapply ( 1:25 , function (x) {  NoFselect [x] = list ( list ( list ( list (DNAm = Allcv_test [[x]] [["DNAm"]] )   ) ,
                                                                  list ( list ( mRNA = Allcv_test [[x]] [["mRNA"]]) )   ,
                                                                  list ( list ( mRNA = Allcv_test [[x]] [["mRNA"]]) ,
                                                                         list ( DNAm = Allcv_test [[x]] [["DNAm"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_test [[x]] [["miRNA"]]) )   ,
                                                                  list ( list ( miRNA = Allcv_test [[x]] [["miRNA"]]) ,
                                                                         list ( DNAm = Allcv_test [[x]] [["DNAm"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_test [[x]] [["miRNA"]]) ,
                                                                         list ( mRNA = Allcv_test [[x]] [["mRNA"]]) ) ,
                                                                  list ( list ( miRNA = Allcv_test [[x]] [["miRNA"]]) ,
                                                                         list ( mRNA = Allcv_test [[x]] [["mRNA"]] )   , 
                                                                         list ( DNAm = Allcv_test [[x]] [["DNAm"]]) )  )  )   } )


#### Noselection 

train_cvs_dataset = t
test_cvs_dataset = s


save (train_cvs_dataset , test_cvs_dataset , labTrainS , labTestS , file = 'CVALLSNFNOS.RData')


save (train_cvs_dataset , test_cvs_dataset , labTrainS , labTestS , file = 'CVALLSNF.RData')


load('CVALLSNF.RData')


SNFPredRT = list ()
for (CV in 1:25) {
  for (MET in 1:6){
    for (OM in 1:7  ){
      DataTrain = lapply  ( train_cvs_dataset [[CV]] [[MET]] [[OM]]  , function (x) x [[1]] )
      DataTest  = lapply  ( test_cvs_dataset  [[CV]] [[MET]] [[OM]]   , function (x) x [[1]] )
      labTrain =  labTrainS  [[CV]] 
      labTest  =  labTestS   [[CV]]
      SNFPredRT = append (SNFPredRT , 
                          list(c (cv     = CV  ,
                                  method      = names( train_cvs_dataset [[CV]] ) [[MET]]  ,
                                  omics       = paste0 (lapply ( test_cvs_dataset [[CV]] [[MET]] [[OM]]  , names  )  , collapse = "_" )  ,
                                  SNF_Pred (DataTrain , DataTest , labTrain , labTest , params , TrueLabs = TRUE)  ) )  )
    }
  }
}


CV = 4
OM = 1 
MET = 7


DataTrain = lapply  ( train_cvs_dataset [[CV]] [[MET]] [[OM]]  , function (x) x [[1]] )
DataTest  = lapply  ( test_cvs_dataset  [[CV]] [[MET]] [[OM]]   , function (x) x [[1]] )
labTrain =  labTrainS  [[CV]] 
labTest  =  labTestS   [[CV]]


###############################################
############ SNF Pred Tests ######################
###############################################


NameSNF =  list.files(path= "01_Results/SNF_R/", 
                      pattern="False.RDS",
                      recursive = TRUE,
                      all.files=TRUE,
                      full.names=FALSE)


Params = lapply (NameSNF , function (x)  ( gsub (".RDS" , "" , x) %>% strsplit ( . , "_") )[[1]] [2:5] )
SNF_R = lapply (NameSNF , function (x) readRDS(  paste0( "01_Results/SNF_R/" , x ) )  )

for (i in 1:length(SNF_R))
{levels (SNF_R[[i]] 
         omics ) == 'c("mRNA", "DNAm")' ] = 'mRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "DNAm")' ] = 'miRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "mRNA")' ] = 'miRNA_mRNA'
levels (SNF_R[[i]] 
        omics ) == 'c("miRNA", "mRNA", "DNAm")' ] = 'miRNA_mRNA_DNAm'}
SNF_R = lapply  (SNF_R , function (D) { cbind ( D  [ , 1:3] , 
                                                D  [, 4:ncol(D) ] %>% apply (. , 2 ,function(x)  as.double ( as.character( x))) ) %>% mutate( CV = rep (1:5 , each =  7*7*5 ) )
  
})






#NameSNF =  list.files(path= "02_Results/SNF_R/", 
#pattern="False.RDS",

NameSNF =  list.files(path= "01_Results/SNF_R_AUC/", 
                      pattern="False_auc.RDS",
                      recursive = TRUE,
                      all.files=TRUE,
                      full.names=FALSE)


Params = lapply (NameSNF , function (x)  ( gsub (".RDS" , "" , x) %>% strsplit ( . , "_") )[[1]] [2:5] )
#SNF_R = lapply (NameSNF , function (x) readRDS(  paste0( "02_Results/SNF_R/" , x ) )  )
SNF_R = lapply (NameSNF , function (x) readRDS(  paste0( "01_Results/SNF_R_AUC/" , x ) )  )
for (i in 1:length(SNF_R))
{levels (SNF_R[[i]] 
         omics ) == 'c("mRNA", "DNAm")' ] = 'mRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "DNAm")' ] = 'miRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "mRNA")' ] = 'miRNA_mRNA'
levels (SNF_R[[i]] 
        omics ) == 'c("miRNA", "mRNA", "DNAm")' ] = 'miRNA_mRNA_DNAm'}
SNF_R = lapply  (SNF_R , function (D) { cbind ( D  [ , 1:3] , 
                                                D  [, 4:ncol(D) ] %>% apply (. , 2 ,function(x)  as.double ( as.character( x))) ) %>% mutate( CV = rep (1:5 , each =  7*7*5 ) )
  
})



SNF_R_DF = SNF_R   %>% do.call( rbind , .)
SNF_R_DF$ID = rep (1:length(SNF_R) , each = nrow(SNF_R[[1]]) )
Params = Params %>% as.data.frame %>% t %>% as.data.frame %>% dplyr ::slice( rep(1:n() , each = nrow(SNF_R[[1]] )) )                
colnames(Params) = c("neighbours","iters","alpha","truelabs") 
SNF_R_DF = cbind (SNF_R_DF , Params )






#NameSNF =  list.files(path= "02_Results/SNF_R/", 
#pattern="False.RDS",

NameSNF =  list.files(path= "01_Results/SNF_NOS_AUC/", 
                      pattern="False_auc.RDS",
                      recursive = TRUE,
                      all.files=TRUE,
                      full.names=FALSE)


Params = lapply (NameSNF , function (x)  ( gsub (".RDS" , "" , x) %>% strsplit ( . , "_") )[[1]] [2:5] )
#SNF_R = lapply (NameSNF , function (x) readRDS(  paste0( "02_Results/SNF_R/" , x ) )  )
SNF_R = lapply (NameSNF , function (x) readRDS(  paste0( "01_Results/SNF_NOS_AUC/" , x ) )  )
for (i in 1:length(SNF_R))
{levels (SNF_R[[i]] 
         omics ) == 'c("mRNA", "DNAm")' ] = 'mRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "DNAm")' ] = 'miRNA_DNAm'
levels (SNF_R[[i]] omics ) == 'c("miRNA", "mRNA")' ] = 'miRNA_mRNA'
levels (SNF_R[[i]] 
        omics ) == 'c("miRNA", "mRNA", "DNAm")' ] = 'miRNA_mRNA_DNAm'}
SNF_R = lapply  (SNF_R , function (D) { cbind ( D  [ , 1:3] , 
                                                D  [, 4:ncol(D) ] %>% apply (. , 2 ,function(x)  as.double ( as.character( x))) ) %>% mutate( CV = rep (1:5 , each =  7*5 ) )
  
})



SNF_S_DF = SNF_R   %>% do.call( rbind , .)
SNF_S_DF$ID = rep (1:length(SNF_R) , each = nrow(SNF_R[[1]]) )
Params = Params %>% as.data.frame %>% t %>% as.data.frame %>% dplyr ::slice( rep(1:n() , each = nrow(SNF_R[[1]] )) )                
colnames(Params) = c("neighbours","iters","alpha","truelabs") 
SNF_S_DF = cbind (SNF_S_DF , Params )


SNF_S_DF$method = "No Selection"


library(ggplot2)
source ("/home//amazigh.mokhtari/NeuroDev_ADD/R/pub_theme.R")


Train_m_auc = rbind (SNF_S_DF, SNF_R_DF )  %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Train_auc = mean (AucTrain)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( auc = mean (Train_auc)) %>% arrange(-auc)

Test_m_auc = rbind (SNF_S_DF, SNF_R_DF )  %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Test_auc = mean (AucTest)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( auc = mean (Test_auc)) 
Train_m_acc = SNF_R_DF %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Train_acc = mean (train_overallAccuracy)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( acc = mean (Train_acc)) %>% arrange(-acc)
Test_m_acc = SNF_R_DF %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Test_acc = mean (test_overallAccuracy)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( acc = mean (Test_acc)) 


Train_m_auc = SNF_R_DF %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Train_auc = mean (AucTrain)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( auc = mean (Train_auc)) %>% arrange(-auc)


Test_m_auc = SNF_R_DF %>% group_by (method , omics , CV , neighbours,iters,alpha) %>% 
  summarize(Test_auc = mean (AucTest)) %>% group_by (method , omics  , neighbours,iters,alpha)   %>% summarize( auc = mean (Test_auc)) 

