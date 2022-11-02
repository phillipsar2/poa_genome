### File for analyses genetic diversity in Poa population panel  ###

# Genome ----
ref = "data/genome/poa_10262021/Poa_pratensis_v1.hap1.fasta"
contig_list = "data/genome/poa_10262021/hap1.contig.list"
fai = "data/genome/poa_10262021/Poa_pratensis_v1.hap1.fasta.fai"
#sample_map = "data/processed/poa.sample_map"

# Process bams ----
sort_in = "data/interm/mapped_bam/{sample}.mapped.bam"
sort_out = "data/interm/sorted_bam/{sample}.sorted.bam"
add_rg = "data/interm/addrg/{sample}.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.dedup.bam"



# Nucleotide diversity --------

## Boulder population
#bamlist = "data/interm/mark_dups/boulder_bamlist.txt"
#saf = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.1dp4.saf.gz"
#prefix = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.1dp4"
#ind = 5
#sfs = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.1dp4.sfs"
#stats = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.1dp4.thetas.idx.pestPG"

## Across populations
bamlist = "data/interm/mark_dups/crosspops_bamlist.txt"
saf = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.1dp4.saf.gz"
prefix = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.1dp4"
ind = 3
sfs = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.1dp4.sfs"
stats = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.1dp4.thetas.idx.pestPG"

# Fst ----

## bamlists for each population (rule saf_per_pop). Choose one at a time.

bamlist_pop  = "data/interm/mark_dups/boulder_bamlist.txt"
ind_pop = 5
saf_pop = "data/saf/boulder/pratensis.boulder.1dp4.saf.gz"
prefix_pop = "data/saf/boulder/pratensis.boulder.1dp4"

#bamlist_pop = "data/interm/mark_dups/argyle_bamlist.txt"
#ind_pop = 1
#saf_pop = "data/saf/argyle/pratensis.argyle.1dp4.saf.gz"
#prefix_pop = "data/saf/argyle/pratensis.argyle.1dp4"

#bamlist_pop = "data/interm/mark_dups/tolstoi_bamlist.txt"
#ind_pop = 1
#saf_pop = "data/saf/tolstoi/pratensis.tolstoi.1dp4.saf.gz"
#prefix_pop = "data/saf/tolstoi/pratensis.tolstoi.1dp4"

# Files for extracting the consensus marker sequences from the Poa population panel ----

## Reference marker sequence for a gene (see Supplementary Materials for list of sequences)
#gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta"

## Process bams
#sort_in = "data/interm/mapped_bam/{sample}.{gene}.mapped.bam"
#sort_out = "data/interm/sorted_bam/{sample}.{gene}.sorted.bam"
#add_rg = "data/interm/addrg/{sample}.{gene}.rg.bam"
#mark_dups = "data/interm/mark_dups/{sample}.{gene}.dedup.bam"

## Call SNPS
#haplo = "data/vcf/{gene}/{sample}.{gene}.vcf"


