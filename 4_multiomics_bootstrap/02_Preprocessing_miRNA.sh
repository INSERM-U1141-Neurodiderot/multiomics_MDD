#!/bin/bash
#MSUB -r 02_Preprocessing_miRNA		# Requested job name
#MSUB -s MDD				# Related study name
#MSUB -o logs/02_Preprocessing_miRNA_${EXPERIMENT_NAME}_.%I.out	# Job standard output (%I is the job ID)
#MSUB -e logs/02_Preprocessing_miRNA_${EXPERIMENT_NAME}_.%I.err		# Job error output (%I is the job ID)
#MSUB -E "-c 24 --mem 220G -t 1440"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

bootstrap_samples_file=$(find . -name "bootstrap_samples_*$EXPERIMENT_NAME.RDS")

module load rstudio

# Start icarust simulation
Rscript 02_Preprocessing_miRNA.R ${bootstrap_samples_file} 24
