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

saveRDS(cv_fold , file = "results/1_CrossValidation/cv_fold.RDS") 
