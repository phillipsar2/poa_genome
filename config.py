### Calling SNPs for whole genome ###

# Genome
ref = "data/genome/poa_01062021/Ppratensis.normalized.nonull.fasta"
contig_list = "data/genome/poa_01062021/Ppratensis.contig.list"
sample_map = "data/processed/poa.sample_map"

# Process bams
sort_out = "data/interm/sorted_bam/{sample}.sorted.bam"
add_rg = "data/interm/addrg/{sample}.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.dedup.bam"



### Consensus sequence building ###

# Mark duplicates
#mark_in = "data/interm/addrg/{sample}.TLF.ref.rg.bam"
#mark_out = "data/interm/mark_dups/{sample}.TLF.ref.dedup.bam"

# Pulling out the chloroplast pseudogenomes
#sort_out = "data/sorted_bam/{sample}.TLF.ref.sorted.bam"
# Aligning a gene to the psuedo ref
#gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta"
#sort_in = "data/interm/mapped_bam/trnTtrnLtrnF/{sample}.matk.mapped.bam"


