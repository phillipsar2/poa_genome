# Genome

ref = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta"
contig_list = "data/gene/trnTtrnLtrnF/trnTtrnLtrnF.contig.list"


# Map reads with bwa mem
bwa_r1 = "data/raw/sequences/{sample}_1.fq.gz",
bwa_r2 = "data/raw/sequences/{sample}_2.fq.gz"

# Mark duplicates
mark_in = "data/interm/addrg/{sample}.TLF.ref.rg.bam"
mark_out = "data/interm/mark_dups/{sample}.TLF.ref.dedup.bam"

# Pulling out the chloroplast pseudogenomes
sort_out = "data/sorted_bam/{sample}.TLF.ref.sorted.bam"
# Aligning a gene to the psuedo ref
gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta"
#sort_in = "data/interm/mapped_bam/trnTtrnLtrnF/{sample}.matk.mapped.bam"



# Multiple sequence alignment
gene_fasta = "data/gene/matK/ALL.matk.fasta"
multi_aln = "data/multi_aln/matk_mafft.fasta"

