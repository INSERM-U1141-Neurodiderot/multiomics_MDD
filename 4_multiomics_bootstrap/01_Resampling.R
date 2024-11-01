library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

cohort          = args[1]
experiment_name = args[2]
B               = args[3]

### input: covariables
covariates = readRDS(file = paste0("data/", cohort, "/covariates.RDS"))
pd_mdd     = readRDS(file = "data/pooled/pd_mdd.RDS")
if (cohort == "female"){
  CovDNAm = c("Array" ,  "Sex" , "BMI.bin")
} else {
  CovDNAm = c("Array" ,  "Sex" , "BMI.bin" , "Age_bin")
}


bootstrap_samples = list()

for (sample in 1:B ) {
  print(sample)
  check_sample = F
  while(!check_sample){
    if (cohort == "pooled"){
      cur_sample = data.frame(x = rownames(covariates), group = interaction(covariates$GROUP, covariates$SEX)) %>% 
        group_by(group)                                                           %>% 
        sample_n(size = n(), replace = T)                                         %>% 
        ungroup 
    }else{
      cur_sample = data.frame(x = rownames(covariates), group = covariates$GROUP) %>% 
        group_by(group)                                                           %>% 
        sample_n(size = n(), replace = T)                                         %>% 
        ungroup 
    }
    tmp          = min(apply(pd_mdd[unique(cur_sample$x), CovDNAm], 2, function(x) min(table(x))))
    check_sample = (tmp > 1)
    print(check_sample)
    if (check_sample){
      if (sample > 1){
        check_duplicated = sapply(bootstrap_samples, function(x){
          count = 0
          for (i in intersect(cur_sample$x, x)){
            tmp = min(sum(cur_sample$x == i), sum(x == i))
            count = count + tmp
          }
          return(count/length(x))
        })
        check_duplicated = max(check_duplicated)
        print(check_duplicated)
        if (check_duplicated < 1){
          bootstrap_samples = append(bootstrap_samples, list(cur_sample$x))
        } else {
          check_sample = FALSE
        }
      } else {
        bootstrap_samples = append(bootstrap_samples, list(cur_sample$x))
      }
    }
  }
}

saveRDS(bootstrap_samples , file = paste0("results/01_Resampling/", cohort, "/bootstrap_samples_", experiment_name, ".RDS"))
