# Rules for analysis of genetic diversity and population structure in Poa population panel

## Alinment ----
# (1) Align pair-end reads to the reference genome
rule bwa_map:
    input:
        ref = config.ref,
        r1 = "data/raw/sequences/{sample}_1.fq.gz",
        r2 = "data/raw/sequences/{sample}_2.fq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t 4 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

# (2) Sort raw bam file
rule samtools_sort:
    input:
        config.sort_in
    output:
        temp(config.sort_out),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

# (3) Add read groups
rule add_rg:
    input:
        config.sort_out
    output:
        bam = temp(touch(config.add_rg))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}",
        sample = "{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk AddOrReplaceReadGroups \
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
        -RGPL=illumina \
        -RGPU=unit1 \
        -RGSM={params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")


# (4) Mark duplicate reads
rule mark_dups:
    input:
        config.add_rg
    output:
        bam = config.mark_dups,
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")

# (5) Evaluate quality metrics with qualimap
rule bamqc:
    input:
        config.mark_dups
    output:
        "reports/bamqc/{sample}_stats/qualimapReport.html"
    params:
        dir = "reports/bamqc/{sample}_stats"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 100000 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=64G")

## Nucleotide diversity ----
# (1) Randomly grab a single read at all sites
rule pi_single:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/pratensis_bamlist.txt",
    output:
        "data/angsd_pi/pratensis.1dp4.MQ20.{run}.ibs.gz"
    params:
        prefix = "data/angsd_pi/pratensis.1dp4.MQ20.{run}"
    run:
        shell("angsd \
        -GL 1 -P 15 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 \
        -setMinDepthInd 1 -setMaxDepthInd 4 \
        -minInd 7 \
        -bam {input.bamlist} \
        -doCounts 1 \
        -doMaf 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}")

# (2) Calculate nucleotide diversity in custom R script

### PCA --------

# Run PCA with single read sampling approach
# -doIBS 1 prints a randomly sampled read from each individual at each position
# -doCov 1 prints out the covariance matrix which can be used for a PCA
# -minMaf 0.05 excludes sites with a minor allele freq less the 0.05
# -doMajorMinor 4 specifies the Major allele to be the reference allele specified by -ref
rule PCA_single:
    input:
        ref = config.ref,
#        bamlist = "data/interm/mark_dups/all_poa_bamlist.txt",
        bamlist = "data/interm/mark_dups/pratensis_bamlist.txt"
    output:
#        "data/pca/all.poa.1dp4.maf05.ibs.gz"
        "data/pca/pratensis.1dp4.maf05.ibs.gz"
    params:
#        prefix = "data/pca/all.poa.1dp4.maf05"
        prefix = "data/pca/pratensis.1dp4.maf05"
    run:
        shell("angsd \
        -GL 1 -P 15 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -setMinDepthInd 1 -setMaxDepthInd 4 \
        -minInd 7 \
        -bam {input.bamlist} \
        -doCounts 1 \
        -doMaf 1 \
        -minMaf 0.05 \
        -ref {input.ref} \
        -doIBS 1 \
        -SNP_pval 1e-4 \
        -doCov 1 \
        -out {params.prefix}")
# make sure to swap minInd back to 8
