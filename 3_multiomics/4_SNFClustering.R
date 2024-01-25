library(caret)
library(SNFtool)
library(pROC)

data_jive_train = readRDS(file = "results/3_FeaturesSelection/data_jive_train.RDS")
data_jive_test = readRDS(file = "results/3_FeaturesSelection/data_jive_test.RDS")
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

cv_all_train = list (miRNA = data_jive_train$miRNA,
                     mRNA  = data_jive_train$mRNA,
                     DNAm  = data_jive_train$DNAm) 
cv_all_test = list (miRNA = data_jive_test$miRNA,
                     mRNA  = data_jive_test$mRNA,
                     DNAm  = data_jive_test$DNAm) 

#####################################################################################

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
  overall   =  ConfMTest$overall
  byClass   =  ConfMTest$byClass
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
  
  # predict labels based on SNF labs of train 
  newLabel =     groupPredict( train =   DataTrain , 
                               test  =   DataTest  ,
                               groups = TrainLab  ,  
                               K = params$K,
                               alpha = params$alpha,
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
  
  predLab =newTestLabel
  TrueLabs=TrueLabTest
  name="test_" 
  
  AriTrain = adjustedRandIndex (TrueLabTrain , TrainLab  )
  AriTest  = adjustedRandIndex (TrueLabTest  , newTestLabel  )
  names(AriTrain) = "AriTrain"
  names(AriTest)  = "AriTest"
  
  #AUC
  AucTrain = as.numeric(roc(as.numeric(TrueLabTrain), as.numeric(TrainLab))$auc)
  AucTest = as.numeric(roc(as.numeric(TrueLabTest), as.numeric(newTestLabel))$auc)
  names(AucTrain) = 'AucTrain'
  names(AucTest) = 'AucTest'
  
  return ( c( ConfMTest , ConfMTrain , AriTrain , AriTest, AucTrain, AucTest)  )
}

labTrain = cov_pooled [cv_all_train$miRNA %>% rownames , "GROUP" ] %>% as.character
labTest = cov_pooled [cv_all_test$miRNA %>% rownames , "GROUP" ] %>% as.character

SNFPredRT = data.frame()
for (t in seq(10, 50, 10)){ #neighbors
  for (K in seq(10, 60, 10)){ #iterations
    for (alpha in seq(0.3, 0.8, 0.1)){ #alpha
      params = list(K = K, alpha = alpha, t = t, Nclust = 2)
      nomParam = paste('K', K, '_alpha', alpha, '_t', t, sep = '')
      print(nomParam)
      for (i in 1:3){
        omics = combn(c('mRNA', 'miRNA', 'DNAm'), i)
        for (om in 1:ncol(omics)){
          DataTrain = list()
          DataTest = list()
          for (o in omics[, om]){
            DataTrain[[o]] = cv_all_train[[o]]
            DataTest[[o]] = cv_all_test[[o]]
          }
          nomOm = paste(omics[, om], collapse = '_')
          temp = data.frame(omics = nomOm, K = params$K, alpha = params$alpha, t = params$t)
          temp = cbind(temp, as.data.frame(t(SNF_Pred (DataTrain, DataTest,
                                                             labTrain, labTest,
                                                             params, TrueLabs = T)))
          )
          if(dim(SNFPredRT)[1] == 0){
            SNFPredRT = temp
          }else{
            SNFPredRT = rbind(SNFPredRT, temp)
          }
          cat(paste(nomOm, SNFPredRT[[nomParam]][[nomOm]]['AriTest'], '\n'))
        }
      }
    }
  }
}

saveRDS(SNFPredRT, 'results/4_SNFClustering/SNFPredRT.RDS')
