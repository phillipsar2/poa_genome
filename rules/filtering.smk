# Extract SNPs

rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/vcf/ETS/{sample}.vcf"
    output:
        "data/raw/vcf_bpres/{sample}.raw.ETS.snps.vcf"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}")


# Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225


# Apply the base GATK filter on the gvcf 

rule hard_filter:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/{sample}.raw.ETS.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{sample}.ETS.filtered.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QD < 2.0\" --filter-name \"QD2\" \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
        -filter \"FS > 60.0\" --filter-name \"FS60\" \
        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
        -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
        -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
        -O {output}")

#######

