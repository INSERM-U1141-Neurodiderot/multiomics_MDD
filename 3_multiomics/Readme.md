# Multi omics analyses
Case-Control multi omics analysis for DNAm, mRNA and miRNA.

## 1 - Cross-validation

A repeated 5-fold CV (5 repetitions) with stratification was used to prevent overfitting. For each re-sampling to be representative of the original dataset, stratification was systematically undertaken based on MDD status, and the Array (DNAm experiments) and sex (pooled cohort only) covariates. For each train set (4/5 of samples), pre-processing, feature selection and clustering were performed as described above. For the test set (remaining 1/5), variables selected on the train set were extracted, pre-processed by applying the different models fitted on the train set and label propagation was used to transfer inferred clusters from the train to the test. For comparison, the same CV procedure was applied to features corresponding to differential analyses results (estimated on each train set separately) and to all features without selection.

- Code: [1_CrossValidation.r](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/1_CrossValidation.r)
- Input: [data/cov_pooled.rds](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data/cov_pooled.rds)
- Output: _cv_fold.rds_, _cv_fold_male.rds_ and _cv_fold_female.rds_ in [results/1_CrossValidation](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/1_CrossValidation)

## 2 - Pre-processing

Pre-processing was applied on each CV folds. In this repository, we provide example for the first folds of pooled, male and female.

- Code: [2_PreProcessing.r](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/2_PreProcessing.r)
- General input: _cov_pooled.rds_, _cv_fold.rds_, _cv_fold_female.rds_ and _cv_fold_male.rds_ in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)

#### mRNA
mRNA pre-processing used DESeq2 library.

- Input: [data/data_mRNA.rds](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data/data_mRNA.rds)
- Output: _cv_mRNA_corr.RDS_, _cv_mRNA_f_corr.RDS_, and _cv_mRNA_m_corr.RDS.rds_ in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

#### miRNA
miRNA pre-processing used DESeq2 library.

- Input: [data/data_miRNA.rds](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data/data_miRNA.rds)
- Output: _cv_miRNA_corr.RDS_, _cv_miRNA_f_corr.RDS_, and _cv_miRNA_m_corr.RDS.rds_ in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

#### DNAm
DNAm pre-processing used neuroComBat v1.0.5 to adjust for covariates, an improved version dedicated to CV procedures: [https://github.com/Jfortin1/ComBatHarmonization/tree/master/R](https://github.com/Jfortin1/ComBatHarmonization/tree/master/R)

- Input: normalised beta-values of probes _myNorm.mdd.RDS_, pd file containes metadata of samples _pd_mdd.RDS_ and leucocyte fractions estimation using Houseman method _LeucocyteFraction.mdd.RDS_ all in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Output: _cv_DNAm_corr.RDS_, _cv_DNAm_f_corr.RDS_, and c_v_DNAm_m_corr.RDS.rds_ in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

In DNAm experiments for the pooled and female cohorts, the Slide covariate was corrected separately in the train and test sets, as it was composed of too many categories for proper representation (in train/test sets) and correction across all CV folds. For similar reasons, in DNAm experiments for the male cohort, both Slide and Array were corrected separately. Due to computational limitations, a pre-filtering was applied to DNAm data: in our own analysis 10% of most variable CpGs selected for downstream analyses. Due to GitHub storage limitation, this frequence was lowed in the code provided 

## 3 - Feature selection with Momix

Building on the [Momix package](https://github.com/cantinilab/momix-notebook) ([Cantini et al. 2021](https://doi.org/10.1038/s41467-020-20430-7)), we encapsulated 6 joint Dimension Reduction (jDR) methods in a CV procedure: RGCCA, JIVE, MCIA, MOFA, intNMF and SciKit-Fusion in order to extract common variance between HDRS matrix scores (considered as a block) and every possible combination of either 1 (mRNA, miRNA, DNAm), 2 (mRNA/miRNA, mRNA/DNAm, miRNA/DNAm) or 3 (mRNA/miRNA/DNAm) omic blocks (7 possible combinations).

For each jDR method and combination of blocks, a factor matrix representing the shared variance across blocks was derived with 10 factors extracted. Among these 10 factors, the one that correlated the most with the MDD status (in absolute value) was kept. Then vectors of weights used to compute this best factor were retrieved. These weights represent the contribution of each feature to the factor to which they belong. For each modality, only the top 10%-ranking elements were retained based on their weigths (in absolute value) of contribution to the best factor. The code is provided for MCIA, MOFA (use python), JIVE, RGCCA and IntNMF jDR methods. The obtained outputs were made available for JIVE using the pre-pocessed data of the first fold of the female cohort as computed in the previous step.

- Code: [3_FeatureSelection.R](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/3_FeatureSelection.R)
- Input: covariables ([data/cov_pooled.RDS](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data/cov_pooled.RDS), and pre-pocessed data of the first fold of the female dataset (_cv_DNAm_f_corr.rds_, _cv_miRNA_f_corr.rds_, _cv_mRNA_f_corr.rds_ in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing))
- Output: _data_jive_test.RDS_ and _data_jive_train.RDS_ in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection)

## 4 - Multiomic clustering with SNF
To evaluate the ability of the selected features to estimate meaningful clusters of subjects in regard to MDD status, we used the integrative clustering technique, Similarity Network Fusion (SNF, [Wang et al. 2014](http://www.nature.com/nmeth/journal/v11/n3/full/nmeth.2810.html)). The code is presented for SNF clustering on the first fold of the female dataset reduced by JIVE method.

- Code: [4_SNFClustering.R](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/4_SNFClustering.R)
- Input: _data_jive_test.RDS_ and _data_jive_train.RDS_ in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection))
- Output: [SNFPredRT.RDS](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/4_SNFClustering/SNFPredRT.RDS)

Several sets of SNF parameters were evaluated through CV
- neighbors ∈ {10, 20, . . . , 50}
- iters ∈ {10, 20, . . . , 60}
- alpha ∈ {0.3, 0.4, . . . , 0.8})

NB: R version 4.3.2 (2023-10-31 ucrt)

