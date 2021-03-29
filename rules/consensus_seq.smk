#### Generate consensus sequences for genes across Poa pop panel ####

## Align ##

# Align a single gene to the chloroplast genome consensus sequences (rpoB-trnC;trnL-trnF)
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


# or #

# Align individuals to a single plastid gene using bwa mem (ITS;ETS;matK)
rule bwa_pull:
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

## Process bams ## (see other rule)

## Calling ##
rule haplotype_caller_gene:
    input:
        ref = config.gene, 
        bam = config.mark_dups
    output:
        config.haplo
#    params:
#        regions = config.contig_list
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output} \
        --reference {input.ref} \
        --G StandardAnnotation \
        -G AS_StandardAnnotation")
#        -L {params.regions}")
#        -ERC BP_RESOLUTION")

## Extract SNPs ##

## Hard filter ## (see other rule)

## Generate consensus sequence ##

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
#rule pseudo_ref:
#    input:
#        ref = config.ref,
#        vcf = "data/processed/filtered_snps_bpres/{sample}.ETS.filtered.snps.vcf"
#    output:
#        "data/processed/pseudo_ref/{sample}.pseudo.fasta"
#        "data/gene/ETS/{sample}.consensus.ETS.fasta"
#    run:
#        shell("gatk IndexFeatureFile -I {input.vcf}")
#        shell("gatk FastaAlternateReferenceMaker \
#        -R {input.ref} \
#        -O {output} \
#       -V {input.vcf}")


# or #

# Pull the gene sequences out of the pseudo-chloroplast genome
#rule get_gene:
#    input:
#        "data/processed/pseudo_ref/{sample}.pseudo.fasta"
#    output:
#        "data/processed/gene/trnTtrnLtrnF/{sample}.trnTtrnLtrnF.fasta"
#    params:
#        pos = 17886,
#        len = 1116
#    rule:
#        shell("~/toolsfordayz/bioawk/bioawk -c fastx ''{{ print $seq }}'' {input} | \
#        awk ''{{print substr($1,17889,1116) }}'' >> {output}")
