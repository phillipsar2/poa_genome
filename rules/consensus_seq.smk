#rule filt_indels:
#    input:
#        vcf = "data/gvcf/{sample}.g.vcf.gz",
#        ref = config.ref
#    output:
#        "data/processed/filt_indels/{sample}.norm.flt-indels.bcf"
#    run:
#        shell("bcftools norm -f {input.ref} {input.vcf} -Ob | \
#        bcftools filter --IndelGap 5 - -Ob -o {output}")
        
# Convert vcfs to iupac format using an older vcf tool
# Input files must be bgzipped and have a tabix index

rule bcf_to_iupac:
    input:
#        bcf = "data/processed/filt_indels/{sample}.norm.flt-indels.bcf"
#        vcf = "data/processed/filtered_snps_bpres/{sample}.filtered.snps.vcf"
        vcf = "data/processed/filtered_snps_bpres/{sample}.filtered.snps.vcf.gz"
    output:
#        gz = "data/processed/filtered_snps_bpres/{sample}.vcf.gz",
        txt = "data/processed/iupac/{sample}.txt",
        head = "data/processed/iupac/{sample}.header"
    run:
#        shell("bcftools view {input.bcf} -v snps -Oz -o {output.vcf}")
        shell("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n' {input.vcf} > {output.txt}")
        shell("bcftools query -f 'chr\tpos\tref\talt[\t%SAMPLE]\n' {input.vcf} > {output.head}")
#        shell("bgzip {input}")
#        shell("tabix -p vcf {output.gz}")
#        shell("vcftools vcf-to-tab {input} > {output.txt}")


# Generate pseudo reference chloroplast genomes with gatk4 FastaAlternateReferenceMaker
# If there are multiple variants that start at a site, it chooses one of them randomly
# This tol only works for SNPs and simple indels
rule pseudo_ref:
    input:
        ref = config.ref,
        vcf = "data/processed/filtered_snps_bpres/{sample}.ETS.filtered.snps.vcf"
    output:
#        "data/processed/pseudo_ref/{sample}.pseudo.fasta"
        "data/gene/ETS/{sample}.consensus.ETS.fasta"
    run:
#        shell("gatk IndexFeatureFile -I {input.vcf}")
        shell("gatk FastaAlternateReferenceMaker \
        -R {input.ref} \
        -O {output} \
        -V {input.vcf}")
