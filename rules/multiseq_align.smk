# Pull the gene sequences out of the pseudo-chloroplast genome
rule get_gene:
    input:
        "data/processed/pseudo_ref/{sample}.pseudo.fasta"
    output:
        "data/processed/gene/trnTtrnLtrnF/{sample}.trnTtrnLtrnF.fasta"
    params:
        pos = 17886,
        len = 1116
    rule:
        shell("~/toolsfordayz/bioawk/bioawk -c fastx ''{{ print $seq }}'' {input} | \
        awk ''{{print substr($1,17889,1116) }}'' >> {output}")


# Aligning the gene sequences using MAFFT
rule multi_aln:
    input:
        config.gene_fasta
    output:
        config.multi_aln
    run:
        shell("mafft \
        --localpair \
        --maxiterate 1000 \
        {input} > {output}")

# Align gene sequences using Clustal Omega
rule clustal_aln:
    input:
        config.gene_fasta
    output:
    run:
        shell("clustalo -i {input} 
