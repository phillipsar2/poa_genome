# Align pair-end reads to the reference genome
rule bwa_map:
    input:
        ref = config.ref,
        r1 = "data/raw/sequences/{sample}_1.fq.gz",
        r2 = "data/raw/sequences/{sample}_2.fq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.ITS.mapped.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"


# Align a single gene to the consensus sequences
rule bwa_single:
    input: 
        gene = config.gene,
        ref = "data/processed/pseudo_ref/{sample}.pseudo.fasta"
    output:
        "data/interm/mapped_bam/{sample}.rpoB-trn-C.mapped.bam"
    run:
        shell("bwa index {input.ref}")
        shell("samtools faidx {input.ref}")
        shell("bwa mem {input.ref} {input.gene} | samtools view -Sb > {output}")

# Takes the input file and stores a sorted version in a different directory.
#rule samtools_sort:
#    input:
#        config.sort_in
#        "data/interm/mapped_bam/{sample}.mapped.bam",
#    output:
#        temp("data/sorted_bam/{sample}.sorted.bam"),
#    params:
#        tmp = "/scratch/aphillip/sort_bam/{sample}"
#    run:
#        shell("mkdir -p {params.tmp}")
#        shell("samtools sort -T {params.tmp} {input} > {output}")
#        shell("rm -rf {params.tmp}")

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
