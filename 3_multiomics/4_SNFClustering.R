library(caret)
library(SNFtool)
library(pROC)
library(mclust)

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

##############################################################################

data_RGCCA_train = readRDS(file = "results/3_FeaturesSelection/data_RGCCA_train_filter_10.RDS")
data_RGCCA_test = readRDS(file = "results/3_FeaturesSelection/data_RGCCA_test_filter_10.RDS")
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

final_results_RGCCA_SNF = list()
for (i in 1:length(data_RGCCA_train)){
  cv_all_train = list (miRNA = data_RGCCA_train[[i]]$miRNA,
                       mRNA  = data_RGCCA_train[[i]]$mRNA,
                       DNAm  = data_RGCCA_train[[i]]$DNAm) 
  cv_all_test = list (miRNA = data_RGCCA_test[[i]]$miRNA,
                      mRNA  = data_RGCCA_test[[i]]$mRNA,
                      DNAm  = data_RGCCA_test[[i]]$DNAm) 
  
  labTrain = cov_pooled [cv_all_train$miRNA %>% rownames , "GROUP" ] %>% as.character
  labTest = cov_pooled [cv_all_test$miRNA %>% rownames , "GROUP" ] %>% as.character
  
  SNFPredRT = data.frame()
  for (t in seq(10, 50, 10)){ #neighbors
    for (K in seq(10, 60, 10)){ #iterations
      for (alpha in seq(0.3, 0.8, 0.1)){ #alpha
        params = list(K = K, alpha = alpha, t = t, Nclust = 2)
        nomParam = paste('K', K, '_alpha', alpha, '_t', t, sep = '')
        print(nomParam)
        nomOm = paste(c('mRNA', 'miRNA', 'DNAm'), collapse = '_')
        temp = data.frame(omics = nomOm, K = params$K, alpha = params$alpha, t = params$t)
        temp = cbind(temp, as.data.frame(t(SNF_Pred (cv_all_train, cv_all_test,
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
  final_results_RGCCA_SNF[[i]] = SNFPredRT
}

saveRDS(final_results_RGCCA_SNF, 'results/4_SNFClustering/SNFPredRT.RDS')

# Display boxplot of AUC for each CV split in the test sets for the optimal set of parameters
res     = sapply(final_results_RGCCA_SNF, function(x) x$AucTest)
res_all = apply(res, 1 , mean)

boxplot(res[which.max(res_all), ], main = "Pooled data - 25 CV split - JDR Method", ylab = "AUC - Test")


stripchart(res[which.max(res_all), ],              # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 4,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 
