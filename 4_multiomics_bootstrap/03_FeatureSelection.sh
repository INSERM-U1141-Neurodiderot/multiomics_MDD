#!/bin/bash
#MSUB -r 03_FeatureSelection            # Requested job name
#MSUB -s MDD                            # Related study name
#MSUB -o logs/03_FeatureSelection_${EXPERIMENT_NAME}_jDR_${1}_omics_${2}.%I.out  # Job standard output (%I is the job ID)
#MSUB -e logs/03_FeatureSelection_${EXPERIMENT_NAME}_jDR_${1}_omics_${2}.%I.err  # Job error output (%I is the job ID)
#MSUB -E "-c 24 --mem 200G -t 1440"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

module load rstudio

# Start icarust simulation
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 1
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 2
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 3
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 4
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 5
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 6
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 7
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 8
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 9
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 10
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 11
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 12
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 13
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 14
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 15
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 16
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 17
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 18
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 19
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 20
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 21
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 22
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 23
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 24
Rscript 03_FeatureSelection.R ${EXPERIMENT_NAME} 24 $1 $2 reference
