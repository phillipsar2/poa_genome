rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = "data/interm/mark_dups/{sample}.dedup.bam"
    output:
        outdir = "data/vcf/{sample}.vcf"
    params:
        regions = "data/genome/poa.contig.list"
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output.outdir} \
        --reference {input.ref} \
#        --G StandardAnnotation \
#        -G AS_StandardAnnotation \
        -L {params.regions}")
#        -ERC BP_RESOLUTION")

