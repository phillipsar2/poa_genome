rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = config.mark_out
    output:
        outdir = "data/vcf/ETS/{sample}.vcf"
    params:
        regions = config.contig_list
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output.outdir} \
        --reference {input.ref} \
#        --G StandardAnnotation \
#        -G AS_StandardAnnotation \
        -L {params.regions}")
#        -ERC BP_RESOLUTION")

