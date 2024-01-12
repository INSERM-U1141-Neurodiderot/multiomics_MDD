R version 4.3.2 (2023-10-31 ucrt)

# Multi omics analyses
Case-Control multi omics analysis for DNAm, mRNA and miRNA.

We encapsulated 6 joint Dimension Reduction (jDR) methods in a CV procedure: RGCCA, JIVE, MCIA, MOFA2, intNMF and SciKit-Fusion in order to extract common variance between HDRS matrix scores (considered as a block) and every possible combination of either 1 (mRNA, miRNA, DNAm), 2 (mRNA/miRNA, mRNA/DNAm, miRNA/DNAm) or 3 (mRNA/miRNA/DNAm) omic blocks.

Snakefile take the omics dataset (in 01_Datasets folder) to compute factorized omics file for evry jDR (in 03_Results folder).

##_Extract_Results.r take factorised omics create buy the snakefile to extract components and features (in 03_Results folder).

01_Datasets contains Male and Female omics datasets
02_Scripts contains 6 jDR scripts
03_Results contains factorized omics for every jDR, with components and features
