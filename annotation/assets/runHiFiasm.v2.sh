#!/bin/bash
if [ "$#" -ne 3 ] ; then
echo "please provide:"
echo -e "\t\t(1) HiFi data (all SMRT cells merged) in fasta format"
echo -e "\t\t(2) Genome name (genus-species) or common name - as single word"
echo -e "\t\t(3) Genome size in bases"
echo "";
echo "runHiFiasm.sh" ;
echo "";
exit 0;
fi
fasta=$1
genome=$2
size=$3
# run hifi
#hifiasm -o ${genome} -t ${SLURM_JOB_CPUS_PER_NODE} -l3 ${fasta}
# convert gfa to gasta
#for ctg in *_ctg.gfa; do
#   awk '$1=="S" {print ">ctg_"NR"\n"$3}' $ctg | fold > ${ctg%.*}.fasta;
#done
hap1=$(find . -name "*hap1.p_ctg.fasta")
hap2=$(find . -name "*hap2.p_ctg.fasta")
primary=$(find . -name "*p_ctg.fasta")
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate ngmlr
# map raw reads back to genome
for asm in $hap1 $hap2 $primary; do
new_Assemblathon.pl -csv -genome_size $size ${asm} > ${asm%.*}.stats
minimap2 -ax map-hifi -t ${SLURM_JOB_CPUS_PER_NODE} ${asm} ${fasta} > ${asm%.*}-hifiasm-ctg_rawreads-mapped.sam
samtools faidx ${asm}
sam=${asm%.*}-hifiasm-ctg_rawreads-mapped.sam
samtools view --threads ${SLURM_JOB_CPUS_PER_NODE} -ht ${asm} -b -o ${sam%.*}.bam ${sam}
samtools sort -o ${sam%.*}_sorted.bam -T ${sam%.*}_temp --threads ${SLURM_JOB_CPUS_PER_NODE} ${sam%.*}.bam
samtools coverage ${sam%.*}_sorted.bam > ${sam%.*}_coverage-stats.txt
ml bbmap
pileup.sh in=${sam} ref=${asm} out=${asm%.*}_counts.txt secondary=f -Xmx100g
awk '$3 > 50000 && $2>=2 {print $1}' ${asm%.*}_counts.txt |grep -v "^#" > ${asm%.*}_cov.ge5x_len.ge50kb.ids
seqtk subseq ${asm}  ${asm%.*}_cov.ge5x_len.ge50kb.ids >  ${asm%.*}_filtered.fasta
done
# clean up
mkdir -p hifiasm-files stats-files
mv *.gfa *.bed *.bin hifiasm-files/
mv *.info *.csv stats-files/
