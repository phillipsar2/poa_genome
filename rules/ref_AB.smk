# calculating allele balance for the reference genome
rule map_pacbio:
    input:
        ref = config.ref,
        pacbio = "data/raw/pacbio/Poa_pratensis_falcon-ec_canu-trimmed_reads.fasta.gz"
    output:
        bam = temp("data/interm/mapped_bam/pacbio.mapped.bam")
    run:
#        shell("module load minimap2/2.16r922")
        shell("minimap2 -ax map-pb {input.ref} {input.pacbio} | \
        samtools view -Sb > {output.bam}")

# sort bam 
# Takes the input file and stores a sorted version in a different directory.
rule sort_pacbioBam:
    input:
        "data/interm/mapped_bam/pacbio.mapped.bam"
    output:
        "data/interm/sorted_bam/pacbio.sorted.bam",
    params:
        tmp = "/scratch/aphillip/sort_bam/pacbio"    
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

# call variable sites in the reference
#rule haplotype_caller:
#    input:
#        ref = config.ref, 
#        bam = "data/interm/sorted_bam/pacbio.sorted.bam"
#    output:
#        outdir = "data/vcf/pacbio.vcf"
#    params:
#        regions = config.contig_list
#    run:
#        shell("gatk HaplotypeCaller \
#        --input {input.bam} \
#        --output {output.outdir} \
#        --reference {input.ref} \
#        --G StandardAnnotation \
#        -G AS_StandardAnnotation \
#        -L {params.regions}")

# call variable sites but with mpileup
rule sam_mpileup:
    input:
        ref = config.ref,
        bam = "data/interm/sorted_bam/pacbio.sorted.bam"
    output:
        vcf = "data/vcf/pacbio.{chrom}.vcf"
    params:
        chr = "{chrom}"
    run:
        shell("module load bcftools")
        # default only sites with max 250 reads considered at each positin, this is way above the max coverage
        # -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
        shell("bcftools mpileup -r {params.chr} -Ou -f {input.ref} {input.bam} \
        --annotate FORMAT/AD,FORMAT/DP | \
        bcftools call -mv -Oz -o {output.vcf}")
