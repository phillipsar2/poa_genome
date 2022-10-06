#!/bin/bash
dir=$1
cd $dir
mkdir -p 1c_purge-haplotigs
cwd=$(pwd)
cpus=${SLURM_JOB_CPUS_PER_NODE}
genome=$(find $(pwd) -name "${dir}.bp.p_ctg.fasta")
bam=$(find $(pwd) -name "${dir}.bp.p_ctg-hifiasm-ctg_rawreads-mapped_sorted.bam")
cd 1c_purge-haplotigs
ln -s $genome
ln -s $bam
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate purge_haplotigs
purge_haplotigs hist -b ${bam} -g ${genome} -t ${cpus}
