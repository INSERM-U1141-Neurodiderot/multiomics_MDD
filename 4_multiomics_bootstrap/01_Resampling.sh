#!/bin/bash
add_date_time=$(date '+%d_%m_%Y__%H_%M_%S')
experiment_name=${2}_cohort_${1}_nbsamples_${3}_${add_date_time}

export EXPERIMENT_NAME=$experiment_name

#MSUB -r 01_Resamping		# Requested job name
#MSUB -s MDD			# Related study name
#MSUB -o logs/01_Resamping_${EXPERIMENT_NAME}_.%I.out	# Job standard output (%I is the job ID)
#MSUB -e logs/01_Resamping_${EXPERIMENT_NAME}_.%I.err		# Job error output (%I is the job ID)
#MSUB -E "-c 1 --mem 1G -t 1"     # Inti node requested; Number of CPUs requested; Requested memory; Time limit

module load rstudio

# Start icarust simulation
Rscript 01_Resampling.R $1 ${EXPERIMENT_NAME} $3
