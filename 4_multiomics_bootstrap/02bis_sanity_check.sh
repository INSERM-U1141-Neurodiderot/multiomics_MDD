#!/bin/bash
#MSUB -r 02bis_sanity_check.R            # Requested job name
#MSUB -s MDD                            # Related study name
#MSUB -o logs/02bis_sanity_check_${EXPERIMENT_NAME}_.%I.out  # Job standard output (%I is the job ID)
#MSUB -e logs/02bis_sanity_check_${EXPERIMENT_NAME}_.%I.err  # Job error output (%I is the job ID)
#MSUB -E "-c 8 --mem 80G -t 1440"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

module load rstudio

# Start icarust simulation
Rscript 02bis_sanity_check.R ${EXPERIMENT_NAME}