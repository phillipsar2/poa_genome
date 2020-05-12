# Genome

ref = "data/gene/ITS/GQ324542.1.fasta"

addrg_in = "data/mergensort/ITS/{sample}.merge.sorted.bam"

mark_in = "data/interm/addrg/{sample}.ITS.rg.bam"
mark_out = "data/interm/addrg/{sample}.ITS.dedup.bam"

# Pulling out the chloroplast pseudogenomes
sort_out = "data/sorted_bam/{sample}.ITS.sorted.bam"
# Aligning a gene to the psuedo ref
gene = "data/gene/ITS/GQ324542.1.fasta"
#sort_in = "data/interm/mapped_bam/trnTtrnLtrnF/{sample}.matk.mapped.bam"

# Multiple sequence alignment
gene_fasta = "data/gene/matK/ALL.matk.fasta"
multi_aln = "data/multi_aln/matk_mafft.fasta"

# Pull out reads that align to a gene
bwa_pull = "data/interm/mapped_bam/ITS/{sample}.mapped.bam"
