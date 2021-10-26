### File for analyses genetic diversity in Poa population panel  ###

# Genome
ref = "data/genome/poa_10262021/Poa_pratensis_v1.fasta.gz"
#contig_list = "data/genome/poa_01062021/Ppratensis.contig.list"
#sample_map = "data/processed/poa.sample_map"

# Process bams
sort_in = "data/interm/mapped_bam/{sample}.mapped.bam"
sort_out = "data/interm/sorted_bam/{sample}.sorted.bam"
add_rg = "data/interm/addrg/{sample}.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.dedup.bam"

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

### Files for extracting the consensus marker sequences from the Poa population panel ###

# Reference marker sequence for a gene (see Supplementary Materials for list of sequences)
gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta"

# Process bams
sort_in = "data/interm/mapped_bam/{sample}.{gene}.mapped.bam"
sort_out = "data/interm/sorted_bam/{sample}.{gene}.sorted.bam"
add_rg = "data/interm/addrg/{sample}.{gene}.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.{gene}.dedup.bam"

# Call SNPS
haplo = "data/vcf/{gene}/{sample}.{gene}.vcf"
