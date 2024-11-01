#!/bin/bash
#MSUB -r 02_Preprocessing_DNAm		# Requested job name
#MSUB -s MDD				# Related study name
#MSUB -o logs/02_Preprocessing_DNAm_${EXPERIMENT_NAME}_.%I.out	# Job standard output (%I is the job ID)
#MSUB -e logs/02_Preprocessing_DNAm_${EXPERIMENT_NAME}_.%I.err		# Job error output (%I is the job ID)
#MSUB -E "-c 30 --mem 360G -t 180"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

bootstrap_samples_file=$(find . -name "bootstrap_samples_*$EXPERIMENT_NAME.RDS")

module load rstudio

# Start icarust simulation
Rscript 02_Preprocessing_DNAm.R ${bootstrap_samples_file} 30 $1
