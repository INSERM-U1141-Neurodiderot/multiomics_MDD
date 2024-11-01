#!/bin/bash
#MSUB -r 04_ComputeIntervals            # Requested job name
#MSUB -s MDD                            # Related study name
#MSUB -o logs/04_ComputeIntervals_${EXPERIMENT_NAME}_jDR_${1}_omics_${2}.%I.out  # Job standard output (%I is the job ID)
#MSUB -e logs/04_ComputeIntervals_${EXPERIMENT_NAME}_jDR_${1}_omics_${2}.%I.err  # Job error output (%I is the job ID)
#MSUB -E "-c 12 --mem 100G -t 180"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

module load rstudio

# Start icarust simulation
Rscript 04_ComputeIntervals.R ${EXPERIMENT_NAME} $1 $2 $3

