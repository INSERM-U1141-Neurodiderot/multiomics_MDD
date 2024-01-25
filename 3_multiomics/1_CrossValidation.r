library(dplyr)
library(stringr)

### input: covariables
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

### create 5 split of 5-fold
cov_female = cov_pooled %>% dplyr::filter (SEX == "F" )
cov_male = cov_pooled %>% dplyr::filter (SEX == "M" )

cv_fold =list ()
for (folds in c(1:5) ) {
    cv = cov_pooled %>%
          group_by(GROUP, SEX) %>% ### generates folds identical in proportion in MDD vs Control and Male vs Female to the original dataset 
          sample_frac(1) %>%
          mutate(fold=rep(1:5, length.out=n())) %>%
          ungroup
    for (d in 1:5)
        {cv_fold = append (cv_fold , list(list(train = cv %>% 
                                                   filter (fold != d) %>% 
                                                   .$Samp  , 
                                               test = cv %>% 
                                                   filter (fold == d) 
                                                   %>% .$Samp ) ) ) }
    }

cv_fold_female =list ()
for (folds in c(1:5) ) {
    cv = cov_female %>%
          group_by(GROUP, SEX) %>% ### generates folds identical in proportion in MDD vs Control and Male vs Female to the original dataset 
          sample_frac(1) %>%
          mutate(fold=rep(1:5, length.out=n())) %>%
          ungroup
    for (d in 1:5)
        {cv_fold_female = append (cv_fold_female , list(list(train = cv %>% 
                                                   filter (fold != d) %>% 
                                                   .$Samp  , 
                                               test = cv %>% 
                                                   filter (fold == d) 
                                                   %>% .$Samp ) ) ) }
    }


cv_fold_male =list ()
for (folds in c(1:5) ) {
    cv = cov_male %>%
          group_by(GROUP, SEX) %>% ### generates folds identical in proportion in MDD vs Control and Male vs male to the original dataset 
          sample_frac(1) %>%
          mutate(fold=rep(1:5, length.out=n())) %>%
          ungroup
    for (d in 1:5)
        {cv_fold_male = append (cv_fold_male , list(list(train = cv %>% 
                                                   filter (fold != d) %>% 
                                                   .$Samp  , 
                                               test = cv %>% 
                                                   filter (fold == d) 
                                                   %>% .$Samp ) ) ) }
    }

saveRDS(cv_fold , file = "results/1_CrossValidation/cv_fold.RDS") 
saveRDS(cv_fold_female , file = "results/1_CrossValidation/cv_fold_female.RDS") 
saveRDS(cv_fold_male , file = "results/1_CrossValidation/cv_fold_male.RDS") 


