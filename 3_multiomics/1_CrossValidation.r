library(dplyr)
library(stringr)

### input: covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

cv_fold   = list()
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

train_sex = test_sex = c()
train_mdd = test_mdd = c()
for (i in 1:25){
  print("TRAIN")
  train_sex = c(train_sex, max(abs(table(cov_pooled[cv_fold[[i]]$train, "SEX"])/sum(table(cov_pooled[cv_fold[[i]]$train, "SEX"])) - table(cov_pooled[, "SEX"])/sum(table(cov_pooled[, "SEX"])))))
  train_mdd = c(train_mdd, max(abs(table(cov_pooled[cv_fold[[i]]$train, "GROUP"])/sum(table(cov_pooled[cv_fold[[i]]$train, "GROUP"])) - table(cov_pooled[, "GROUP"])/sum(table(cov_pooled[, "GROUP"])))))
  print("TEST")
  test_sex = c(test_sex, max(abs(table(cov_pooled[cv_fold[[i]]$test, "SEX"])/sum(table(cov_pooled[cv_fold[[i]]$test, "SEX"])) - table(cov_pooled[, "SEX"])/sum(table(cov_pooled[, "SEX"])))))
  test_mdd = c(test_mdd, max(abs(table(cov_pooled[cv_fold[[i]]$test, "GROUP"])/sum(table(cov_pooled[cv_fold[[i]]$test, "GROUP"])) - table(cov_pooled[, "GROUP"])/sum(table(cov_pooled[, "GROUP"])))))
}
max(train_sex)
max(train_mdd)
max(test_sex)
max(test_mdd)

print()
print(table(cov_pooled[, "GROUP"])/sum(table(cov_pooled[, "GROUP"])))
unique(sapply(cv_fold, function(x) length(unique(unlist(x)))))
sapply(cv_fold, function(x) sapply(x, length))

saveRDS(cv_fold , file = "results/1_CrossValidation/cv_fold.RDS") 

