#### Filtering vcfs for poa pop panel #####
# select only bialleleic SNPs
rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/vcf/mpileup/all.poa.raw.vcf.gz" 
    output:
        "data/vcf/mpileup/all.poa.snps.vcf"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        --restrict-alleles-to BIALLELIC \
        -O {output}")


# Filtering diagnostics
# Extract variant quality scores
# https://evodify.com/gatk-in-non-model-organism/

rule diagnostics:
    input:
        vcf = "data/vcf/mpileup/all.poa.snps.vcf",
        ref = config.ref
    output:
        "reports/filtering/all.poa.snps.vcf.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -F AB -F MQM -F MQMR -F AF -GF AO \
        -O {output}")

# apply hard filtering following GATK best practices
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

rule hard_filter:
    input:
        ref = config.ref,
        vcf = "data/vcf/mpileup/all.poa.snps.vcf"
    output:
        filt =  "data/processed/filtered_snps/all.poa.filtered.snps.vcf",
        exclude = "data/processed/filtered_snps/all.poa.filtered.nocall.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 40.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
        -O {output.filt}")
        shell("gatk SelectVariants -V {output.filt} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output.exclude}")

# Evaluate depth across samples to set DP filter

rule depth:
    input:
        vcf = "data/processed/filtered_snps/all.poa.filtered.nocall.snps.vcf",
        ref = config.ref
    output:
        dp = "reports/filtering/all.poa.depth.filtered.nocall.table"
    run:
#        shell("tabix -p vcf {input.vcf}")
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -GF DP \
        -O {output.dp}")

# Fitlter by depth of each individual then filter on maximum nocall

rule filter_depth:
    input:
        vcf = "data/processed/filtered_snps/all.poa.filtered.nocall.snps.vcf",
        ref = config.ref
    output:
        dp = "data/processed/filtered_snps/all.poa.filtered.nocall.2dp20.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 2 || DP > 20 \" \
        -G-filter-name \"DP_2-20\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}")

#### Filtering to building consensus sequences ####

# Extract SNPs

#rule get_snps:
#    input:
#        ref = config.ref,
#        vcf = "data/vcf/ETS/{sample}.vcf"
#    output:
#        "data/raw/vcf_bpres/{sample}.raw.ETS.snps.vcf"
#    run:
#        shell("gatk SelectVariants \
#        -R {input.ref} \
#        -V {input.vcf} \
#        -select-type SNP \
#        -O {output}")


# Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225


# Apply the base GATK filter on the gvcf 

#rule hard_filter:
#    input:
#        ref = config.ref,
#        vcf = "data/raw/vcf_bpres/{sample}.raw.ETS.snps.vcf"
#    output:
#        "data/processed/filtered_snps_bpres/{sample}.ETS.filtered.snps.vcf"
#    run:
#        shell("gatk VariantFiltration \
#        -V {input.vcf} \
#        -filter \"QD < 2.0\" --filter-name \"QD2\" \
#        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
#        -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
#        -filter \"FS > 60.0\" --filter-name \"FS60\" \
#        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
#        -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
#        -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
#        -O {output}")

#######

