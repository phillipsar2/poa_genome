#!/bin/bash
#SBATCH -N 1
#SBATCH -A triffid
#SBATCH -n 36
#SBATCH --mem=350GB
#SBATCH -t 21-00:00:00
#SBATCH -J bionano_0
#SBATCH -o bionano_0.o%j
#SBATCH -e bionano_0.e%j
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module purge
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/gcc/7.3.0
#module use /work/GIF/software/modules
./runBionano.sh poa-v1.fasta EXP_REFINEFINAL1.cmap Reimported_AJ608_Andropogon_gerardii_Filter300_AJ608_Andropogon_gerardii_Filter300_RawMolecules.filtered.bnx
scontrol show job $SLURM_JOB_ID
