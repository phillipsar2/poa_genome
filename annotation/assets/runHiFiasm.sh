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
hifiasm -o ${genome} -t ${SLURM_JOB_CPUS_PER_NODE} -l2 ${fasta}
# convert gfa to gasta
for ctg in *_ctg.gfa; do
   awk '$1=="S" {print ">"$2"\n"$3}' $ctg | fold > ${ctg%.*}.fasta;
done
# combine haplotypes
alt=$(find . -name "*a_ctg.fasta")
pri=$(find . -name "*p_ctg.fasta")
if [[ -f ${pri} && -f ${alt} ]]; then
cat ${alt} ${pri} >> ${genome}_diploid.fasta
else
hap=$(find . -name "*_ctg.fasta")
cp $hap ${genome}_diploid.fasta
fi
# map raw reads back to genome
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate ngmlr
minimap2 -ax map-hifi -t ${SLURM_JOB_CPUS_PER_NODE} ${genome}_diploid.fasta ${fasta} > ${genome}_diploid-hifiasm-ctg_rawreads-mapped.sam
samtools faidx ${genome}_diploid.fasta
sam=${genome}_diploid-hifiasm-ctg_rawreads-mapped.sam
samtools view --threads ${cpus} -ht ${genome}_diploid.fasta.fai -b -o ${sam%.*}.bam ${sam}
samtools sort -o ${sam%.*}_sorted.bam -T ${sam%.*}_temp --threads ${cpus} ${sam%.*}.bam
samtools coverage ${sam%.*}_sorted.bam > ${sam%.*}_coverage-stats.txt
conda deactivate
ml bbmap
pileup.sh in=${genome}_diploid-hifiasm-ctg_rawreads-mapped.sam ref=${genome}_diploid.fasta out=counts.txt secondary=f -Xmx100g
awk '$3 > 50000 && $2>=2 {print $1}' counts.txt |grep -v "^#" > coverage-ge5x_and_len-ge50kb_contigs.ids
seqtk subseq ${genome}_diploid.fasta coverage-ge5x_and_len-ge50kb_contigs.ids > ${genome}_diploid-filtered.fasta
# assembly stats
if [[ -f ${pri} && -f ${alt} ]]; then
new_Assemblathon.pl -csv -genome_size $size ${pri}
new_Assemblathon.pl -csv -genome_size $size ${alt}
else
new_Assemblathon.pl -csv -genome_size $size ${hap}
fi
for fasta in ${genome}_diploid.fasta ${genome}_diploid-filtered.fasta; do
new_Assemblathon.pl -csv -genome_size ${size} ${fasta}
done
for csv in *.csv; do
datamash transpose --field-separator="," < $csv | head -n 30 | awk -F"," '$1 !~ /^Percentage/ {print $2}' > ${csv%.*}.info
done
# clean up
mkdir -p hifiasm-files stats-files
mv *.gfa *.bed *.bin hifiasm-files/
mv *.info *.csv stats-files/
