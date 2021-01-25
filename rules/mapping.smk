# Align pair-end reads to the reference genome
rule bwa_map:
    input:
        ref = config.ref,
        r1 = "data/raw/sequences/{sample}_1.fq.gz",
        r2 = "data/raw/sequences/{sample}_2.fq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"


# Align a single gene to the chloroplast genome consensus sequences
#rule bwa_single:
#    input: 
#        gene = config.gene,
#        ref = "data/processed/pseudo_ref/{sample}.pseudo.fasta"
#    output:
#        "data/interm/mapped_bam/{sample}.rpoB-trn-C.mapped.bam"
#    run:
#        shell("bwa index {input.ref}")
#        shell("samtools faidx {input.ref}")
#        shell("bwa mem {input.ref} {input.gene} | samtools view -Sb > {output}")

# Align a single gene to the poa reference genome
#rule bwa_ref:
#    input:
#        gene = config.gene,
#        ref = "data/genome/{sample}.fasta"
#    output:
#        temp("data/interm/mapped_bam/{sample}.TLF.ref.mapped.bam")
#    run:
#        shell("bwa mem {input.ref} {input.gene} | samtools view -Sb > {output}")

# Align induviduals to a single plastid gene using bwa mem
# Not in use
#rule bwa_pull:
#    input:
#        ref = config.gene,
#        r1 = "data/raw/sequences/{sample}_1.fq.gz",
#        r2 = "data/raw/sequences/{sample}_2.fq.gz"
#    output:
#        config.bwa_pull
#        "data/interm/mapped_bam/{sample}.mapped.bam"
#    log:
#        "logs/bwa_mem/{sample}.log"
#    shell:
#        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
#        "samtools view -Sb > {output}) 2> {log}"
