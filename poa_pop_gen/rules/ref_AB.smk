# Calculating allele balance for the reference genome


# (1) Align error corrected PacBio reads to the super scaffolds
rule map_pacbio:
    input:
        ref = config.ref,
        pacbio = "data/raw/pacbio/Poa_pratensis_falcon-ec_canu-trimmed_reads.fasta.gz"
    output:
        bam = temp("data/interm/mapped_bam/pacbio.mapped.bam")
    run:
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

# Quality metrics with qualimap
rule bamqc:
    input:
        "data/interm/sorted_bam/pacbio.sorted.bam"
    output:
        "reports/bamqc/pacbio/qualimapReport.html"
    params:
        dir = "reports/bamqc/pacbio"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 100000 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=64G")


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
        # default only sites with max 250 reads considered at each positin, this is way above the max coverage
        # -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
        shell("bcftools mpileup -r {params.chr} -Ou -f {input.ref} {input.bam} \
        --annotate FORMAT/AD,FORMAT/DP | \
        bcftools call -mv -Oz -o {output.vcf}")

# merge the reference scaffold vcfs for ease of use
rule merge_vcfs:
    input:
        "data/vcf/pacbio_variant_files.list"
    output:
        "data/vcf/pacbio.all_scaff.vcf.gz"
    run:
        shell("bcftools concat -f {input} -O z -o {output}")

# Extract GT and DP for calculating AB
rule diag_ab:
    input:
        vcf = "data/vcf/pacbio.all_scaff.vcf.gz",
        ref = config.ref
    output:
        "reports/filtering/pacbio.AB.table"
    run:
        shell("tabix -p vcf {input.vcf}")
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -GF GT -GF AD \
        -O {output}")

# calculate allele balance (AB) at every site
rule calc_AB:
    input:
        "reports/filtering/pacbio.AB.table"
    output:
        "reports/filtering/pacbio.AB.estimate.txt"
#    params:
#        "1"
#    conda:
#        "env/r4.yaml"
    shell:
        "Rscript scripts/allelebalance_filter.R"
