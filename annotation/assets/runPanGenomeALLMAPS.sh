#!/bin/bash
if [ "$#" -lt 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) geneome fasta file"
echo -e "\t\t(2) short name"
echo -e "\t\t(3) list of scafs to exclude (txt file, 1 id per line) (optional)"
echo "";
echo "./runPanGenomeALLMAPS.sh <genome.fasta> <name> [exlude.txt]" ;
echo "";
exit 0;
fi
genome="$1"
name="$2"
cpus="${SLURM_JOB_CPUS_PER_NODE:-36}"
ml purge
module load samtools
module load bedtools2
module load hisat2
remove=$(find $(pwd) -name "remove.txt")
# check if the genome is softlinked
if [ -L "$genome" ]; then
    echo "the genome fasta file is softlinked. the container will not be able to see the directory it is located"
    echo "attempting to make a local copy"
    cp $(realpath ${genome}) $(basename ${genome%.*}).copy.fasta
    genome=$(basename ${genome%.*}).copy.fasta 
else
    echo "local genome file."
fi

# check if index exists
check=$(find $(pwd)/index -name "*.ht2*" | head -n 1)
if [ -f "$check" ]; then
    echo "index for this genome exists, skipping building the index.."
else
    echo "no index files found, running HiSat2 indexing now.."
    mkdir -p index
    hisat2-build -p $cpus $genome index/${name}
fi
index="$(pwd)/index/${name}"
markers="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/pb_anchors.fasta"
info="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/pb_anchors.txt"
# map markers
hisat2 -p ${cpus} --mp 1,1 --no-softclip -f -x ${index} -U ${markers}  1> ${name}_mapped.sam 2> mapping_stats.txt
# sort and covert to bam
samtools view -b -o ${name}_mapped.bam ${name}_mapped.sam
samtools sort -o ${name}_sorted.bam ${name}_mapped.bam
samtools view -h ${name}_sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > unique_mapped.bam
mv unique_mapped.bam ${name}_sorted.bam
bedtools bamtobed -i ${name}_sorted.bam |  awk '{print $4"\t"$1"\t"$2}' > ${name}_part1.txt
# generate info
awk 'BEGIN{OFS=FS="\t"}{print $2,$1,$3}' ${info} > ${name}_part2.txt
echo "Scaffold ID,scaffold position,LG,genetic position" > ${name}_pg-anchor.csv
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' ${name}_part2.txt ${name}_part1.txt | sed 's/ //g' | cut -f 2- |sed 's/\t/,/g' | sort | uniq >> ${name}_pg-anchor.csv
# cleanup markers
ml purge
ml bioawk
marker=${name}_pg-anchor.csv
rm  pg-markers-mapped.csv &> /dev/null
bioawk -c fastx '{print $name"\t"length($seq)}' ${genome} > scaf-sizes.txt
grep -v "^Scaffold ID," ${marker} | cut -f 1 -d "," | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1 -n > ${name}-marker-density.txt
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' scaf-sizes.txt ${name}-marker-density.txt |awk '{print $0"\t"($2/$3)*10000}' > ${name}-marker-table.txt
awk '$3>100000 && $2>20' ${name}-marker-table.txt | cut -f 1 > scafs-for-agp.txt
# exclude the scafs that are in remove file (if it exists)
if [ -f "$3" ]; then
    echo "remove file exists, removing those ids from map"
    grep -Fvw -f $3 scafs-for-agp.txt > temp
    mv temp scafs-for-agp.txt
fi
# subset the markers to make it run faster
while read line;  do
   awk -v x=$line -F"," '$1==x'  ${name}_pg-anchor.csv | sort | uniq | sort -R | head -n 100;
done<scafs-for-agp.txt >> pangenome.csv

# optional filtering steps:
cut -f 1,3 pangenome.csv -d "," | sed 's/,/\t/g' | sort |uniq -c | awk '{print $2"\t"$3"\t"$1}' |sort -k1,1 -k3,3rn  | awk '!seen[$1]++ {print "ID:"$1",chr"$2}' > needed.txt
awk -F"," '{print "ID:"$1",chr"$3"\t"$0}' pangenome.csv > temp.csv
grep -Fw -f needed.txt temp.csv |cut -f 2 > pangenome_filtered.csv
mv pangenome.csv pangenome_orig.csv
mv pangenome_filtered.csv pangenome.csv
# run ALLMAPS
ml purge
ml singularity
scmd="singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/jcvi.sif"
$scmd python3 -m jcvi.assembly.allmaps merge pangenome.csv -o pangenome.bed
$scmd python3 -m jcvi.assembly.allmaps path --mincount=3 --format=png pangenome.bed ${genome}
