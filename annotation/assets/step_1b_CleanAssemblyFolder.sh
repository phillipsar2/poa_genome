#!/bin/bash
dir=$1
cd $dir
mkdir -p 0_data/{assembly,rawdata} 1_hifiasm 2_genome-stats/{coverage-stats,busco,assembly-metrics}
mv *.fasta.gz ./0_data/rawdata
mv *.stats 2_genome-stats/assembly-metrics/
mv hifi.sub runHiFiasm.sh nova* 1_hifiasm/
mv hifiasm-files/* 1_hifiasm/
rmdir hifiasm-files
mv stats-files/*.csv 2_genome-stats/assembly-metrics/
rmdir stats-files
mv *.txt 2_genome-stats/coverage-stats/
rm *ctg_rawreads-mapped.bam *ctg_rawreads-mapped.sam *ge50kb.ids *_filtered.fasta *.fai
mv *p_ctg.fasta 0_data/assembly
mv *.bam 1_hifiasm/
