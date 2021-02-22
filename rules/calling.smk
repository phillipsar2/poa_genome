######## bcftools SNP calling pipeline #############

# calling SNPs for the population panel then calculating quality and depth statistics for filtering
rule sam_mpileup:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/bamlist.txt"
    output:
        # outputs a gvcf but bcftools treats like vcf (?)
        bcf = "data/bcf/all.poa.raw.bcf",
    run:
        shell("module load bcftools")
        # default only sites with max 250 reads considered at each positin, this is way above the max coverage
        # -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
        shell("bcftools mpileup -Ou -f {input.ref} -b {input.bamlist} \
        --annotate FORMAT/AD,FORMAT/DP | \
        bcftools call -mv -Ob -o {output.bcf}")

rule process_bcf:
    input: 
        bcf = "data/bcf/all.poa.raw.bcf",
        ref = config.ref
    output:
        gz = "data/vcf/mpileup/all.poa.raw.vcf.gz",
        stats = "reports/mpileup/all.poa.raw.vcf.stats"
    params:
        vcf = "data/vcf/mpileup/all.poa.raw.vcf",
    run:
        shell("bcftools view {input.bcf} > {params.vcf}")
        shell("bgzip {params.vcf}")
        shell("gatk IndexFeatureFile -I {output.gz}")
#        shell("bcftools tabix -p vcf {output.gz}")
        shell("bcftools stats -F {input.ref} -s - {input.bcf} > {output.stats}")
        shell("plot-vcfstats -p reports/filtering/plots/ {output.stats}")

######### GATK SNP calling pipeline ##############
# Calling SNPs for the consensus sequences

rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = config.mark_dups
    output:
        outdir = "data/vcf/{gene}.{sample}.vcf"
    params:
        regions = config.contig_list
    run:
        shell("gatk HaplotypeCaller \
        --input {input.bam} \
        --output {output.outdir} \
        --reference {input.ref} \
        --G StandardAnnotation \
        -G AS_StandardAnnotation \
        -L {params.regions} \
        -ERC BP_RESOLUTION")

# Scatter reference into intervals using SplitIntervals. Set intervals in Snakefile.
# https://gatk.broadinstitute.org/hc/en-us/articles/360036348372-SplitIntervals
# ``--scatter-count n` splits reference into n intervals
rule split_intervals:
    input:
        ref = config.ref
    output:
        int = expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
    params:
        regions = config.contig_list,
        dir = "data/processed/scattered_intervals"
    run:
        shell("gatk SplitIntervals -R {input.ref} -L {params.regions} --scatter-count 200 -O {params.dir}")

# Combine GVCFs with GenomicsDBImport
# sample_map file is samplename with path/to/gvcf in the following line
# scattered interval file created by SpitIntevals

rule combine_gvcfs:
    input:
        gvcfs = expand("data/gvcf/{sample}.g.vcf", sample = SAMPLE),
        region = "data/processed/scattered_intervals/{interval}-scattered.interval_list",
        map = config.sample_map
    output:
        directory("data/interm/combined_database_bpres/{interval}")
    params:
        tmp = "/scratch/aphillip/genomicsdbimport/{interval}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk --java-options \"-Xmx90g -Xms90g\" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output} \
        --batch-size 50 \
        --reader-threads 8 \
        --sample-name-map {input.map} \
        --intervals {input.region} --tmp-dir {params.tmp}")
        shell("rm -rf {params.tmp}")


rule joint_geno:
    input:
        dir = directory("data/interm/combined_database_bpres/{interval}"),
        ref = config.ref
    output:
        "data/raw/vcf_bpres/{interval}.raw.vcf"
    params:
        db = "gendb://data/interm/combined_database_bpres/{interval}",
        region = "data/processed/scattered_intervals/{interval}-scattered.interval_list",
    run:
        shell("gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {params.db} \
        -L {params.region} \
        -new-qual \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --include-non-variant-sites \
        -O {output}")

