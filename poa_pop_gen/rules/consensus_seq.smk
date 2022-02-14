#### Generate consensus sequences for genes across Poa pop panel ####

# (1)  Align
rule bwa_map:
    input:
        ref = config.gene,
        r1 = "data/raw/sequences/{sample}_1.fq.gz",
        r2 = "data/raw/sequences/{sample}_2.fq.gz"
    output:
        config.sort_in
    log:
        "logs/bwa_pull/{sample}.log"
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

## (2-5) Process bams ## (see process_bam.smk)

# (6)  Call snps
rule haplotype_caller_gene:
    input:
        ref = config.gene, 
        bam = config.mark_dups
    output:
        "data/vcf/{gene}/{sample}.{gene}.vcf"
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output} \
        --reference {input.ref} \
        --G StandardAnnotation")

# (7) Extract SNPs
rule get_snps:
    input:
        ref = config.gene,
        vcf = "data/vcf/{gene}/{sample}.{gene}.vcf"
    output:
        "data/vcf/{gene}/{sample}.{gene}.snps.vcf"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}")

# (8)  Hard filter
rule hard_filter:
    input:
        ref = config.gene,
        vcf = "data/vcf/{gene}/{sample}.{gene}.snps.vcf"
    output:
        filt =  "data/processed/filtered_snps/{gene}/{sample}.{gene}.filtered.snps.vcf",
        exclude = "data/processed/filtered_snps/{gene}/{sample}.{gene}.filtered.nocall.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 40.0\" --filter-name \"QUAL40\" \
        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
        -O {output.filt}")
        shell("gatk SelectVariants -V {output.filt} --exclude-filtered true  -O {output.exclude}")

# (9)  Generate pseudo reference chloroplast genomes with gatk4 FastaAlternateReferenceMaker
# If there are multiple variants that start at a site, it chooses one of them randomly
# This tol only works for SNPs and simple indels
rule pseudo_ref:
    input:
        ref = config.gene,
        vcf = "data/processed/filtered_snps_bpres/{gene}/{sample}.{gene}.filtered.nocall.snps.vcf"
    output:
        "data/gene/{gene}/{sample}.consensus.{gene}.fasta"
    run:
        shell("gatk IndexFeatureFile -I {input.vcf}")
        shell("gatk FastaAlternateReferenceMaker \
        -R {input.ref} \
        -O {output} \
       -V {input.vcf}")
