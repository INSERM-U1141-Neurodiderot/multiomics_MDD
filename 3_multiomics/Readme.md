R version 4.3.2 (2023-10-31 ucrt) with packages:
- caret (6.0.94)
- ChAMP (2.32.0)
- DESeq2 (1.42.0)
- dplyr (1.1.4)
- IntNMF (1.2.0)
- MOFA2 (1.12.0)
- omicade4 (1.42.0)
- pROC (1.18.5)
- r.jive (2.4)
- reticulate (1.20)
- RGCCA (3.0.3)
- SNFtool (2.3.1)
- stringr (1.5.5)
- sva (3.50.0)

# Multi omics analyses
We implemented a supervised multiomic integration framework designed to identify blood-derived molecular features that best discriminate MDD patients from healthy controls. Our framework was encapsulated in a Cross-Validation (CV) procedure.
In this repository, an example is provided using pooled male and female data together with the JDR method RGCCA. The code for this Case-Control multi omics analysis for DNAm, mRNA and miRNA is presented in several parts, organised as:


## 1 - Cross-validation

A repeated 5-fold CV (5 repetitions) with stratification was used to prevent overfitting. For each re-sampling, we at first wanted to stratify on MDD status, sex but also on the Slide covariate.

However, the number of sub-categories between all these covariates being too high, this stratification was not possible. Not taking into account the Slide covariate when re-sampling was not an option as it would result in categories of Slide with too few observations in the train set to properly correct for the effect of this covariate.

As a result, we decided to still stratify on the MDD status and the sex, and to assign each whole category of Slide either on the train or on the test set. Eventually, we had to correct separately for the Slide effect on the train and on the test set.

For each train set (4/5 of samples), pre-processing, feature selection and clustering were performed. For the test set (remaining 1/5), variables selected on the train set were extracted, pre-processed by applying the different models fitted on the train set (except for the correction of the Slide covariate) and label propagation was used to transfer inferred clusters from the train to the test. 

- Code: [1_CrossValidation.r](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/1_CrossValidation.r)
- Input: _cov_pooled.rds_ in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Generated output placed in [results/1_CrossValidation](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/1_CrossValidation)

## 2 - Pre-processing

Pre-processing was applied on each CV folds (pooled data).

- Code: [2_PreProcessing.r](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/2_PreProcessing.r)
- General input: _cov_pooled.rds_ and _cv_fold.rds_ in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)

#### mRNA
mRNA pre-processing used DESeq2 library.

- Input: _data_mRNA.rds_ in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Generated output placed in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

#### miRNA
miRNA pre-processing used DESeq2 library.

- Input: data_miRNA.rds in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Generated output placed in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

#### DNAm
DNAm pre-processing used _neuroComBat_ v1.0.5 to adjust for covariates, an improved version dedicated to CV procedures: [https://github.com/Jfortin1/ComBatHarmonization/tree/master/R](https://github.com/Jfortin1/ComBatHarmonization/tree/master/R). This part uses the _reticulate_ library, which requires python3.

- Input: pd file containes metadata of samples _pd_mdd.RDS_ and leucocyte fractions estimation using Houseman method _LeucocyteFraction.mdd.RDS_ all in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data). Beta-values of probes _GSE251780_Mokhtari_MatrixSignalIntensities.txt_ availiable at [geo](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251786)
- Generated output placed in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

## 3 - Feature selection with Momix

Building on the [Momix package](https://github.com/cantinilab/momix-notebook) ([Cantini et al. 2021](https://doi.org/10.1038/s41467-020-20430-7)), we encapsulated 6 joint Dimension Reduction (jDR) methods in a CV procedure: RGCCA, JIVE, MCIA, MOFA and intNMF and SciKit-Fusion in order to extract common variance between HDRS matrix scores (considered as a block) and every possible combination of either 1 (mRNA, miRNA, DNAm), 2 (mRNA/miRNA, mRNA/DNAm, miRNA/DNAm) or 3 (mRNA/miRNA/DNAm) omic blocks (7 possible combinations).

For each jDR method and combination of blocks, a factor matrix representing the shared variance across blocks was derived with 10 factors extracted. Among these 10 factors, the one that correlated the most with the MDD status (in absolute value) was kept. Then vectors of weights used to compute this best factor were retrieved. These weights represent the contribution of each feature to the factor to which they belong. For each modality, only the top 10%-ranking elements were retained based on their weigths (in absolute value) of contribution to the best factor. The code is provided for RGCCA. Code is also available for MCIA, MOFA, JIVE, and IntNMF jDR methods, but commented.

- Code: [3_FeatureSelection.R](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/3_FeatureSelection.R)
- Input: covariables (_cov_pooled.RDS_ in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)), and pre-pocessed data of the first fold of the dataset (_cv_DNAm_corr.rds_, _cv_miRNA_corr.rds_, _cv_mRNA_corr.rds_ previously generated in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing))
- Generated output placed in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection)

## 4 - Multiomic clustering with SNF
To evaluate the ability of the selected features to estimate meaningful clusters of subjects in regard to MDD status, we used the integrative clustering technique, Similarity Network Fusion (SNF, [Wang et al. 2014](http://www.nature.com/nmeth/journal/v11/n3/full/nmeth.2810.html)). The code of SNF clustering is presented for pooled dataset reduced by RGCCA method. Several sets of SNF parameters were evaluated through CV:
--- neighbors ∈ {10, 20, . . . , 50}
--- iters ∈ {10, 20, . . . , 60}
--- alpha ∈ {0.3, 0.4, . . . , 0.8}

- Code: [4_SNFClustering.R](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/4_SNFClustering.R)
- Input generated in 'feature selection' steps (_data_test.RDS_ and _data_train.RDS_) and saved in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection)
- Generated output placed in: [results/4_SNFClustering](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/4_SNFClustering)

