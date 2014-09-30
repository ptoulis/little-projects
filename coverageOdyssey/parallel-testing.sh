#!/bin/bash
#SBATCH -J parallel-testing # name for job array
#SBATCH -o all.out #Standard output
#SBATCH -e all.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 00:20:00 #Running time of 20 mins.
#SBATCH --mem-per-cpu 3000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine

# first arg = true mean (unknown)
# second arg = sd (known)
# third arg = no. of samples (N)
# fourth arg = number of reps / node
# fifth argument = job id
Rscript main-odyssey.R 2.5 1.5 100 50 $SLURM_ARRAY_TASK_ID
