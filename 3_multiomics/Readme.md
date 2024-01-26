R version 4.3.2 (2023-10-31 ucrt)

# Multi omics analyses
Case-Control multi omics analysis for DNAm, mRNA and miRNA.

## 1 - Cross-validation

A repeated 5-fold CV (5 repetitions) with stratification was used to prevent overfitting.

For each re-sampling to be representative of the original dataset, stratification was systematically undertaken based on MDD status, and the Slide (DNAm experiments) and sex (pooled cohort only) covariates.

For each train set (4/5 of samples), pre-processing, feature selection and clustering were performed as described above.

For the test set (remaining 1/5), variables selected on the train set were extracted, pre-processed by applying the different models fitted on the train set and label propagation was used to transfer inferred clusters from the train to the test.

For comparison, the same CV procedure was applied to features corresponding to differential analyses results (estimated on each train set separately) and to all features without selection.

- [1_CrossValidation](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/1_CrossValidation.r)
- Input: cov_pooled.rds in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Output: cv_fold.rds, cv_fold_male.rds and cv_fold_female.rds in [results/1_CrossValidation](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/1_CrossValidation)

## 2 - Pre-processing

Pre-processing was applied on each CV folds. In this repository, we provide example for the first folds of pooled, male and female.

- [2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/2_PreProcessing.r)
- General input: cov_pooled.rds, cv_fold.rds, cv_fold_female.rds and cv_fold_male.rds in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)

#### mRNA
mRNA pre-processing used DESeq2 library.

- Input: data_mRNA.rds in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Output: cv_mRNA_corr.RDS, cv_mRNA_f_corr.RDS, and cv_mRNA_m_corr.RDS.rds in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

#### miRNA
miRNA pre-processing used DESeq2 library.

- Input: data_miRNA.rds in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Output: cv_miRNA_corr.RDS, cv_miRNA_f_corr.RDS, and cv_miRNA_m_corr.RDS.rds in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)


#### DNAm
DNAm pre-processing used neuroComBat v1.0.5 to adjust for covariates, an improved version dedicated to CV procedures: https://github.com/Jfortin1/ComBatHarmonization/tree/master/R

In DNAm experiments for the pooled and female cohorts, the Slide covariate was corrected separately in the train and test sets, as it was composed of too many categories for proper representation (in train/test sets) and correction across all CV folds.

For similar reasons, in DNAm experiments for the male cohort, both Slide and Array were corrected separately.

- Input: normalised beta-values of probes (myNorm.mdd.RDS), pd file containes metadata of samples (pd_mdd.RDS) and leucocyte fractions estimation using Houseman method (LeucocyteFraction.mdd.RDS) all in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data)
- Output: cv_DNAm_corr.RDS, cv_DNAm_f_corr.RDS, and cv_DNAm_m_corr.RDS.rds in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

Due to computational limitations, a pre-filtering was applied to DNAm data: in our own analysis 10% of most variable CpGs selected for downstream analyses.
Due to GitHub storage limitation, this frequence was lowed in the code provided 

## 3 - Feature selection with Momix

Building on the Momix benchmark, we encapsulated 6 joint Dimension Reduction (jDR) methods in a CV procedure: RGCCA, JIVE, MCIA, MOFA, intNMF and SciKit-Fusion in order to extract common variance between HDRS matrix scores (considered as a block) and every possible combination of either 1 (mRNA, miRNA, DNAm), 2 (mRNA/miRNA, mRNA/DNAm, miRNA/DNAm) or 3 (mRNA/miRNA/DNAm) omic blocks (7 possible combinations).

For each jDR method and combination of blocks, a factor matrix representing the shared variance across blocks was derived with 10 factors extracted.

Among these 10 factors, the one that correlated most (in absolute value) with MDD status was kept.

Then weight vectors used to compute the best factor were retrieved. They represent the contribution of each feature to the factor to which they belong. Thus, for each modality, only the top 10%-ranking elements of weight vectors (in absolute value) were retained.

- Input are covariables (cov_pooled) in [data](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/data), and pre-pocessed data of female fold 1 computed in the previous step (cv_DNAm_f_corr.rds, cv_miRNA_f_corr.rds, cv_mRNA_f_corr.rds) in [results/2_PreProcessing](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/2_PreProcessing)

- Output in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection)

- Code was provided for MCIA, MOFA (use python), JIVE, RGCCA and IntNMF: [3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/3_FeatureSelection.r). However, only JIVE features were actually computed here.

## 4 - Multiomic clustering with SNF
To evaluate the ability of selected features to estimate meaningful clusters of subjects in regard to MDD status, we used an unsupervised integrative clustering technique, SNF.

Here we provided code for SNF clustering on the first fold female dataset reduced by JIVE method: [4_SNFClustering](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/4_SNFClustering.r), using data from previous step (in [results/3_FeatureSelection](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/3_FeatureSelection)). Results in [results/4_SNFClustering](https://github.com/INSERM-U1141-Neurodiderot/multiomics_MDD/tree/main/3_multiomics/results/4_SNFClustering)

Several sets of SNF parameters were evaluated through CV
- neighbors ∈ {10, 20, . . . , 50}
- iters ∈ {10, 20, . . . , 60}
- alpha ∈ {0.3, 0.4, . . . , 0.8})

