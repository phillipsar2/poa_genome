rule haplotype_caller:
    input:
        ref = config.ref, 
        bam = config.mark_dups
    output:
        outdir = "data/vcf/{sample}.vcf"
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

