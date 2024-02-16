library(dplyr)
library(stringr)

### input: covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

### create 5 split of 5-fold
cov_female = cov_pooled %>% dplyr::filter (SEX == "F" )
cov_male = cov_pooled %>% dplyr::filter (SEX == "M" )

cv_fold = list()
new_GROUP = as.factor(apply(table(interaction(cov_pooled$GROUP, cov_pooled$SEX), cov_pooled$Slide), 2, 
                            function(x) paste0(x, collapse = "_")))
for (folds in c(1:5) ) {
  cv = data.frame(x = unique(cov_pooled$Slide), new_GROUP = new_GROUP) %>%
    group_by(new_GROUP) %>% ### generates folds identical in proportion in MDD vs Control and Male vs Female to the original dataset 
    sample_frac(1) %>%
    mutate(fold=rep(1:5, length.out=n())) %>%
    ungroup
  for (i in 1:5){
    slide_to_keep = cv$x[which(cv$fold != i)]
    cv_fold = append(cv_fold,
                     list(list(train = cov_pooled$Samp[which(cov_pooled$Slide %in% slide_to_keep)],
                               test = cov_pooled$Samp[-which(cov_pooled$Slide %in% slide_to_keep)])))
  }
}


for (i in 1:25){
  print("TRAIN")
  print(table(cov_pooled[cv_fold[[i]]$train, "SEX"])/sum(table(cov_pooled[cv_fold[[i]]$train, "SEX"])))
  print(table(cov_pooled[cv_fold[[i]]$train, "GROUP"])/sum(table(cov_pooled[cv_fold[[i]]$train, "GROUP"])))
  print(unique(table(cov_pooled[cv_fold[[i]]$train, "Slide"])))
  print("TEST")
  print(table(cov_pooled[cv_fold[[i]]$test, "SEX"])/sum(table(cov_pooled[cv_fold[[i]]$test, "SEX"])))
  print(table(cov_pooled[cv_fold[[i]]$test, "GROUP"])/sum(table(cov_pooled[cv_fold[[i]]$test, "GROUP"])))
  print(unique(table(cov_pooled[cv_fold[[i]]$test, "Slide"])))
}


print(table(cov_pooled[, "SEX"])/sum(table(cov_pooled[, "SEX"])))
print(table(cov_pooled[, "GROUP"])/sum(table(cov_pooled[, "GROUP"])))
unique(sapply(cv_fold, function(x) length(unique(unlist(x)))))
sapply(cv_fold, function(x) sapply(x, length))







cv_fold_female = list()
new_GROUP_female = as.factor(apply(table(interaction(cov_female$GROUP, cov_female$SEX), cov_female$Slide), 2,
                                   function(x) paste0(x, collapse = "_")))
for (folds in c(1:5) ) {
  cv = data.frame(x = unique(cov_female$Slide), new_GROUP_female = new_GROUP_female) %>%
    group_by(new_GROUP_female) %>% ### generates folds identical in proportion in MDD vs Control and Male vs Female to the original dataset 
    sample_frac(1) %>%
    mutate(fold=rep(1:5, length.out=n())) %>%
    ungroup
  for (i in 1:5){
    slide_to_keep = cv$x[which(cv$fold != i)]
    cv_fold_female = append(cv_fold_female,
                            list(list(train = cov_female$Samp[which(cov_female$Slide %in% slide_to_keep)],
                                      test = cov_female$Samp[-which(cov_female$Slide %in% slide_to_keep)])))
  }
}



cv_fold_female = list()
new_GROUP_female = as.factor(apply(table(interaction(cov_female$GROUP, cov_female$SEX), cov_female$Slide), 2,
                                   function(x) paste0(x, collapse = "_")))
for (folds in c(1:5) ) {
  cv = data.frame(x = unique(cov_female$Slide), new_GROUP_female = new_GROUP_female) %>%
    group_by(new_GROUP_female) %>% ### generates folds identical in proportion in MDD vs Control and Male vs Female to the original dataset 
    sample_frac(1) %>%
    mutate(fold=rep(1:5, length.out=n())) %>%
    ungroup
  for (i in 1:5){
    slide_to_keep = cv$x[which(cv$fold != i)]
    cv_fold_female = append(cv_fold_female,
                            list(list(train = cov_female$Samp[which(cov_female$Slide %in% slide_to_keep)],
                                      test = cov_female$Samp[-which(cov_female$Slide %in% slide_to_keep)])))
  }
}

for (i in 1:25){
  print("TRAIN")
  print(table(cov_female[cv_fold_female[[i]]$train, "SEX"])/sum(table(cov_female[cv_fold_female[[i]]$train, "SEX"])))
  print(table(cov_female[cv_fold_female[[i]]$train, "GROUP"])/sum(table(cov_female[cv_fold_female[[i]]$train, "GROUP"])))
  print(unique(table(cov_female[cv_fold_female[[i]]$train, "Slide"])))
  print("TEST")
  print(table(cov_female[cv_fold_female[[i]]$test, "SEX"])/sum(table(cov_female[cv_fold_female[[i]]$test, "SEX"])))
  print(table(cov_female[cv_fold_female[[i]]$test, "GROUP"])/sum(table(cov_female[cv_fold_female[[i]]$test, "GROUP"])))
  print(unique(table(cov_female[cv_fold_female[[i]]$test, "Slide"])))
}

print(table(cov_female[, "SEX"])/sum(table(cov_female[, "SEX"])))
print(table(cov_female[, "GROUP"])/sum(table(cov_female[, "GROUP"])))
unique(sapply(cv_fold_female, function(x) length(unique(unlist(x)))))
sapply(cv_fold_female, function(x) sapply(x, length))





cv_fold_male = list()
new_GROUP_male = as.factor(apply(table(interaction(cov_male$GROUP, cov_male$SEX), cov_male$Slide), 2,
                                 function(x) paste0(x, collapse = "_")))
for (folds in c(1:5) ) {
  cv = data.frame(x = unique(cov_male$Slide), new_GROUP_male = new_GROUP_male) %>%
    group_by(new_GROUP_male) %>% ### generates folds identical in proportion in MDD vs Control and Male vs male to the original dataset 
    sample_frac(1) %>%
    mutate(fold=rep(1:5, length.out=n())) %>%
    ungroup
  for (i in 1:5){
    slide_to_keep = cv$x[which(cv$fold != i)]
    cv_fold_male = append(cv_fold_male,
                          list(list(train = cov_male$Samp[which(cov_male$Slide %in% slide_to_keep)],
                                    test = cov_male$Samp[-which(cov_male$Slide %in% slide_to_keep)])))
  }
}

for (i in 1:25){
  print("TRAIN")
  print(table(cov_male[cv_fold_male[[i]]$train, "SEX"])/sum(table(cov_male[cv_fold_male[[i]]$train, "SEX"])))
  print(table(cov_male[cv_fold_male[[i]]$train, "GROUP"])/sum(table(cov_male[cv_fold_male[[i]]$train, "GROUP"])))
  print(unique(table(cov_male[cv_fold_male[[i]]$train, "Slide"])))
  print("TEST")
  print(table(cov_male[cv_fold_male[[i]]$test, "SEX"])/sum(table(cov_male[cv_fold_male[[i]]$test, "SEX"])))
  print(table(cov_male[cv_fold_male[[i]]$test, "GROUP"])/sum(table(cov_male[cv_fold_male[[i]]$test, "GROUP"])))
  print(unique(table(cov_male[cv_fold_male[[i]]$test, "Slide"])))
}

print(table(cov_male[, "SEX"])/sum(table(cov_male[, "SEX"])))
print(table(cov_male[, "GROUP"])/sum(table(cov_male[, "GROUP"])))
unique(sapply(cv_fold_male, function(x) length(unique(unlist(x)))))
sapply(cv_fold_male, function(x) sapply(x, length))

saveRDS(cv_fold , file = "cv_fold.RDS") 
saveRDS(cv_fold_female , file = "cv_fold_female.RDS") 
saveRDS(cv_fold_male , file = "cv_fold_male.RDS") 

