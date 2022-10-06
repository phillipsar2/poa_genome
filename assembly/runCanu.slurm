#!/bin/bash
#SBATCH -N 1
#SBATCH -A triffid
#SBATCH -n 36
#SBATCH -p compute
#SBATCH -t 1-00:00:00
#SBATCH -J falcon
#SBATCH -o falcon.o%j
#SBATCH -e falcon.e%j
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module purge
PATH=$PATH:/work/LAS/mhufford-lab/shared_dir/minconda/20181213/bin
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate denovo_asm
canu -trim-assemble -p pratensis -d canu-20200118 -s canu.cfg -pacbio-corrected consensus.merged.fa.gz
jobid=$(squeue -o "%8i" -u arnstrm |sort -n |tail -n 1)
sleep 1800
sbatch --dependency=afterok:${jobid} canu.sub
scontrol show job $SLURM_JOB_ID
