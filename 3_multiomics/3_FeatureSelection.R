library(dplyr)
library(omicade4)
library(r.jive)
library(RGCCA)
library(MOFA2)
library(IntNMF)

cv_DNAm_f_corr = readRDS("results/2_PreProcessing/cv_DNAm_corr.RDS")
cv_miRNA_f_corr = readRDS("results/2_PreProcessing/cv_miRNA_corr.RDS")
cv_mRNA_f_corr = readRDS("results/2_PreProcessing/cv_mRNA_corr.RDS")
cov_pooled = readRDS(  file = "data/cov_pooled.RDS")

omics_test = list(miRNA = cv_miRNA_f_corr$corrected_miRNA_test,
                  mRNA = cv_mRNA_f_corr$corrected_mRNA_test,
                  DNAm = cv_DNAm_f_corr$corrected_DNAm_test,
                  covariates = cov_pooled[rownames(cv_miRNA_f_corr$corrected_miRNA_test),])
omics_train = list(miRNA = cv_miRNA_f_corr$corrected_miRNA_train,
                  mRNA = cv_mRNA_f_corr$corrected_mRNA_train,
                  DNAm = cv_DNAm_f_corr$corrected_DNAm_train,
                  covariates = cov_pooled[rownames(cv_miRNA_f_corr$corrected_miRNA_train),])

omics_test$covariates$GROUP = as.numeric (as.factor (omics_test$covariates$GROUP))
# omics_test$covariates = omics_test$covariates %>% select(GROUP)
omics_test$covariates = as.matrix (omics_test$covariates)

omics_train$covariates$GROUP = as.numeric (as.factor (omics_train$covariates$GROUP))
# omics_train$covariates = omics_train$covariates %>% select (GROUP)
omics_train$covariates = as.matrix (omics_train$covariates)

rm(cv_DNAm_f_corr)
rm(cv_miRNA_f_corr)
rm(cv_mRNA_f_corr)
rm(cov_pooled)

# Number of Factors
num.factors = 10

################
###   MCIA   ###
################

MCIA_dry = mcia(lapply(omics_train[1:3], t)  , cia.nf = num.factors)

# Extract factors
factors_mcia = as.matrix(MCIA_dry$mcoa$SynVar)
colnames (factors_mcia) = paste0 ( "MCIA_", 1:10)

# Extract Components 
compo_mcia = lapply(  MCIA_dry$coa , function (x) {  df = x$co ; colnames(df) = paste0("MCIA_Comp_" , 1:10) ; df } )

# Extract features 
features_w_mcia = lapply (MCIA_dry$coa , function (x) { df = x$li ; colnames(df) = paste0("MCIA_" , 1:10) ; df  })

################
###   JIVE   ###
################

factorizations_jive = jive( lapply (omics_train[1:3], t) , 
                            rankJ   = num.factors, 
                            rankA   = rep ( num.factors, length( omics_train ) ), 
                            method  = "given", 
                            conv    = "default", 
                            maxiter = 1000, 
                            showProgress = T)
rankJV = factorizations_jive$rankJ
rankIV.v = factorizations_jive$rankA

J  = numeric(0)
ng = 0

factors_jive = list()
features_w_jive = list()

for(j in 1:length(omics_train [1:3])){
  J = rbind(J, factorizations_jive$joint[[j]])
  ng = c(ng,dim(factorizations_jive$joint[[j]])[1])
}
svd.o = svd(J)
jV = svd.o$v %*% diag(svd.o$d)

for(j in 1:length(omics_train[1:3])){
  features_w_jive[[j]]  =  svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]
  rownames(features_w_jive[[j]]) = rownames( t(omics_train[[j]]))
  colnames(features_w_jive[[j]]) = paste0("JIVE_" , 1:num.factors )
}
factors_jive = jV[, 1:rankJV]

colnames(factors_jive) = paste0("JIVE_" , 1:10 )
rownames(factors_jive) = rownames(omics_train$mRNA)
names (features_w_jive) = c("miRNA",  "mRNA" , "DNAm")

get_jiv_omic_comp = function (factorizations_jive , x) {
  svd_s = svd (factorizations_jive$joint[[x]])
  s_v = svd_s$v %*% diag(svd_s$d)
  s = s_v [,1:10]
  colnames (s) = paste0("JIVE_Compo_" , 1:10)
  rownames (s) = omics_train$covariates %>% rownames
  return(s)
}
compo_jive = lapply (1:3 , function (x) {get_jiv_omic_comp (factorizations_jive , x )} ) 
names(compo_jive) = c("miRNA","mRNA","DNAm")

rm(J)
rm(jV)
rm(svd.o)
rm(j)
rm(ng)
rm(rankIV.v)
rm(rankJV)

#################
###   RGCCA   ###
#################
result_sgcca = rgcca(blocks  = lapply (omics_train[1:3], t)  ,
                     ncomp = rep(num.factors, length(lapply (omics_train[1:3], t))),
                     scheme = "centroid",
                     scale = TRUE,
                     init = "svd",
                     bias = TRUE,
                     tol = 1e-08,
                     verbose = T,
                     superblock = T)

features_w_rgcca = result_sgcca$Y[1:3]
features_w_rgcca$miRNA = features_w_rgcca$miRNA[colnames(omics_train$miRNA),]
features_w_rgcca$mRNA = features_w_rgcca$mRNA[colnames(omics_train$mRNA),]
features_w_rgcca$DNAm = features_w_rgcca$DNAm[colnames(omics_train$DNAm),]

features_w_rgcca = lapply (features_w_rgcca , function (x) { colnames(x) = paste0("RGCCA_" , 1:10) ; x  })

factors_RGCCA = result_sgcca$a$superblock
colnames (factors_RGCCA) = paste0 ( "RGCCA_", 1:10)


##################
###   IntNMF   ###
##################
omics_pos<-list()
for(j in 1:3){
  if(min(omics_train[[j]])<0){
    omics_pos[[j]]<-omics_train[[j]]+abs(min(omics_train[[j]]))
  }else{
    omics_pos[[j]]<-omics_train[[j]]
  }
  omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
}

factorizations_intnmf = nmf.mnnals(dat = lapply (omics_pos, as.matrix ) ,  k  = num.factors)

factors_intnmf = factorizations_intnmf$W
colnames(factors_intnmf) = paste0( "IntNMF_" , 1:10 )

features_w_intNMF =list()

for(j in 1:length(omics_train [1:3])){
  features_w_intNMF[[j]] = t(factorizations_intnmf$H[[j]]) 
  # rownames(features_w_intNMF [[j]] ) = rownames(omics[[j]])
  colnames(features_w_intNMF [[j]] ) = paste0( "IntNMF_" , 1:10 )
}
names (features_w_intNMF) = c("miRNA",  "mRNA" , "DNAm")

rm(omics_pos)

################
###   MOFA   ###
################

# Creating MOFA object from a list of matrice
MOFAobject = create_mofa( lapply (omics_train[1:3] , t ) )
DataOptions = get_default_data_options(MOFAobject)
model_opts = get_default_model_options(MOFAobject)
model_opts$num_factors = num.factors 
TrainOptions = get_default_training_options(MOFAobject)

# Prepare the MOFA Object 
MOFAobject =  prepare_mofa( object = MOFAobject,
                            data_options = DataOptions,
                            model_options = model_opts,
                            training_options = TrainOptions )

# Configure python for MOFA2
use_basilisk = FALSE
reticulate::use_python("/opt/miniconda3/bin/python3.7", required=TRUE)
outfile = file.path(getwd(),"model2.hdf5")

MOFAobject_trained = run_mofa(MOFAobject, outfile)

features_w_mofa = get_weights(MOFAobject_trained)
features_w_mofa = lapply (features_w_mofa , function(x) {colnames (x) = paste0 ("MOFA_" , 1:10) ; x } )

factors_mofa = get_factors (MOFAobject_trained) $group1 
colnames (factors_mofa) = paste0("MOFA_" , 1:10)

compo_mofa = MOFAobject_trained@expectations$W
compo_mofa = lapply (compo_mofa , function (x) {colnames(x) = paste0("MOFA_Compo" , 1:10) ; x } )

#########################
###   Best Features   ###
#########################

### get the ranking of the feature based on the most correlated feature to the phenotype

## get best corr to phenotype 
Group = as.double(as.factor(omics_train$covariates[, 'GROUP']))
best_cor = colnames(factors_jive)[which.max(abs(apply(factors_jive, 2, function(x){cor(x, Group)})))]

gen_omic_n = function (features_w , best_fact , o , omics , frac = 0.1) {
  if ( !is.null (best_fact) ) { # get row names of the omic
    feat_n =  features_w[[ o ]]
    feat_n =  feat_n %>% as.data.frame
    feat_n =  feat_n %>% arrange ( - abs( . [ , best_fact ]))
    feat_n =  feat_n %>% top_frac (  frac )
    feat_n =  feat_n %>% row.names
    return( omics [[ o]] [ , feat_n ]  )
  }
}

fracts = list (miRNA = 0.1 , mRNA = 0.1 , DNAm =  0.1)
om = list (miRNA =  "miRNA" , mRNA =  "mRNA" , DNAm = "DNAm")

data_jive_train = mapply(function(x , y) {gen_omic_n (features_w_jive , best_cor , x , omics = omics_train , frac = y )}  , om , fracts  )

data_jive_test = omics_test
data_jive_test$mRNA = data_jive_test$mRNA[, colnames(data_jive_train$mRNA)]
data_jive_test$miRNA = data_jive_test$miRNA[, colnames(data_jive_train$miRNA)]
data_jive_test$DNAm = data_jive_test$DNAm[, colnames(data_jive_train$DNAm)]

saveRDS(data_jive_train, 'results/3_FeaturesSelection/data_jive_train.RDS')
saveRDS(data_jive_test, 'results/3_FeaturesSelection/data_jive_test.RDS')

