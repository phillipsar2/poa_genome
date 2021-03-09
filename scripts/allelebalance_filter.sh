#!/bin/bash -l

#SBATCH --job-name=allelebalance
#SBATCH -D /group/jrigrp10/aphillip/poa
#SBATCH -e /group/jrigrp10/aphillip/poa/slurm_log/sterror_%j.txt
#SBATCH -o /group/jrigrp10/aphillip/poa/slurm_log/stdoutput_%j.txt
#SBATCH -p med2
#SBATCH -t 7-00:00
#SBATCH --mem 16G

cat scripts/allelebalance_filter.R

module load R

srun Rscript scripts/allelebalance_filter.R
