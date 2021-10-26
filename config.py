### Calling SNPs for whole genome ###

# Genome
ref = "data/genome/poa_10262021/Poa_pratensis_v1.fasta.gz"
#contig_list = "data/genome/poa_01062021/Ppratensis.contig.list"
#sample_map = "data/processed/poa.sample_map"

# Process bams
## Process population panel
#sort_in = "data/interm/mapped_bam/{sample}.mapped.bam"
#sort_out = "data/interm/sorted_bam/{sample}.sorted.bam"
#add_rg = "data/interm/addrg/{sample}.rg.bam"
#mark_dups = "data/interm/mark_dups/{sample}.dedup.bam"

# Evaluating allele balance and depth - all Poa
#ab_table = "reports/filtering/all.poa.AB.table"

# Evaluating allele balance and depth - Ppratensis
ab_table = "reports/filtering/pPratensis.AB.table"

## Consensus sequence building --------

# Pulling out the chloroplast pseudogenomes
#sort_out = "data/sorted_bam/{sample}.TLF.ref.sorted.bam"

# bwa_pull - aligning individual's reads to the reference sequence
gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta" # reference
sort_in = "data/interm/mapped_bam/{sample}.trnLtrnF.mapped.bam"

sort_out = "data/interm/sorted_bam/{sample}.trnLtrnF.sorted.bam"
add_rg = "data/interm/addrg/{sample}.trnLtrnF.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.trnLtrnF.dedup.bam"

haplo = "data/vcf/trnLtrnF/{sample}.trnLtrnF.vcf"

# Nucleotide diversity --------

## Boulder population
#bamlist = "data/interm/mark_dups/boulder_bamlist.txt"
#saf = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.2dp6.saf.gz"
#prefix = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.2dp6"
#ind = 5
#sfs = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.2dp6.sfs"
#stats = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.2dp6.thetas.idx.pestPG"

## Across populations
bamlist = "data/interm/mark_dups/crosspops_bamlist.txt"
saf = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.2dp6.saf.gz"
prefix = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.2dp6"
ind = 3
sfs = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.2dp6.sfs"
stats = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.2dp6.thetas.idx.pestPG"


# PCA from snps --------
## convert vcf to beagle
final_vcf = "data/processed/filtered_snps/all.poa.filtered.nocall.2dp20.max0.snps.vcf"
beagle = "data/beagle/all-poa/all.poa.{chrom}.BEAGLE.PL"

## generate beagle file for pca
#bamlist = "data/interm/mark_dups/all_poa_bamlist.txt.txt"
#saf = "data/pca/all.poa.{chrom}.saf.gz"
#prefix = "data/pca/all.poa.{chrom}"
#ind = 8

