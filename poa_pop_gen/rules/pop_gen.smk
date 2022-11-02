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
# modified from script shared by Silas Tittes

## (1) Generating saf file - site allele frequency likelihood
# -GL 1: samtools method for GL, request 15 threads,
# -uniqueOnly 1: unique reads only; -remove_bads: remove reads with flag above 255; only proper pairs; -trim 0: don't trim, -C 50: adjust mapQ for excessive mismatches (as SAMtools)
# -minMapQ 30: minimum mapping qual 30; -minQ 30: minimum quality score 30
# -doCounts: calculate the frequency of different bases (required for depth filter
# -setMinDepthInd 2 -setMaxDepthInd 20: depth of an individual at a site must be 2 <= x <= 20
# -minInd: num of individuals required to have data at a site (we specificy 100% of samples)
# -doSaf 1: generate saf file, calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
        # -doMaf = estimate allele frequency
        # -doMajorMinor = how to decide the major allele, 4: for major alle according to reference states

rule angsd_saf:
    input:
        ref = config.ref,
        bamlist = config.bamlist
    output:
        config.saf
    params:
        ind = {config.ind},
        prefix = config.prefix,
        chrom = "{chrom}"
    run:
        shell("angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 \
        -setMinDepthInd 1 -setMaxDepthInd 4 \
        -minInd {params.ind} \
        -ref {input.ref} -anc {input.ref} \
        -doSaf 1 \
        -r {params.chrom} \
        -bam {input.bamlist} \
        -out {params.prefix}")


## (2) Generate global estimate of SFS (site frequency spectrum) from saf file
# do not need to specify region as each saf is a region
# realSFS generates the maximum likelihood estimate of the SFS
# -fold 1 specifies the folded spectrum as I don't have an ancestral state
### when using a folded SFS, only thetaW, thetaD, and tajimasD will be meaningful in the output of realSFS
# real SFS saf2theta calculates the thetas for each site
rule pop_sfs:
    input:
        ref = config.ref,
        saf = config.saf
    output:
        sfs = config.sfs
    params:
        prefix = config.prefix
    run:
        shell("realSFS {params.prefix}.saf.idx  -P 10 -fold 1 > {params.prefix}.sfs")
        shell("realSFS saf2theta {params.prefix}.saf.idx -sfs {params.prefix}.sfs -outname {params.prefix}")

# (3) calculate thetas (and neutrality tests) in 10k sliding windows
### when using a folded SFS, only thetaW (tW), thetaD (tP), and tajimasD will be meaningful in the output of realSFS
rule pop_pi:
    input:
        sfs = config.sfs
    output:
        stats = config.stats
    params:
        prefix = config.prefix,
        win = 10000
    run:
        shell("thetaStat do_stat {params.prefix}.thetas.idx \
        -win {params.win} \
        -step {params.win}")
       # merge the files
#        shell("cat <(cat *.thetas.idx.pestPG | head -n1) <(cat *.thetas.idx.pestPG | grep -v nSites) > all.boulder.thetas.idx.pestPG")

# (4) merge files
rule merge_thetas:
    output:
        "data/angsd_pi/{site}/{site}.1dp4.thetas.idx.pestPG"
    params:
        site = "{site}"
    shell:
        """        
        cat <(cat data/angsd_pi/{params.site}/*1dp4.thetas.idx.pestPG | head -n1) <(cat data/angsd_pi/{params.site}/*1dp4.thetas.idx.pestPG | grep -v nSites) > {output}
        """
    

# (1b) Randomly grab a single read at all sites
rule pi_single:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/pratensis_bamlist.txt",
    output:
        "data/angsd_pi/pratensis.1dp4.maf05.MQ20.ibs.gz"
    params:
        prefix = "data/angsd_pi/pratensis.1dp4.maf05.MQ20"
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
        -minMaf 0.05 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}")


### PCA --------

# (1) generate GL in beagle format
# doMajorMinor 4: use refence allele as major
# assuming a mean coverage of 2X, 99% of sites should have 6X coverage in a poission distribution; ppois(6, lambda = 2)
rule angsd_beagle:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/all_poa_bamlist.txt"
    output:
        "data/pca/all.poa.{chrom}.1dp4.beagle.gz"
    params:
        prefix = "data/pca/all.poa.{chrom}.1dp4",
        chrom = "{chrom}"
    run:
        shell("angsd -GL 1 -P 15 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 \
        -setMinDepthInd 1 -setMaxDepthInd 4 \
        -minInd 8 \
        -ref {input.ref} -anc {input.ref} \
        -r {params.chrom} \
        -bam {input.bamlist} \
        -out {params.prefix}")

# (2) merge beagle files together
#        shell("cat <(zcat data/pca/*beagle.gz | head -n1) <(zcat data/pca/*beagle.gz | \
#        grep -v marker) > data/pca/all.poa.merged.2dp6.be

# (3a) run PCA

rule angsd_pca:
    input:
        beagle = expand("data/pca/all.poa.{chrom}.1dp4.beagle.gz", chrom = CHROM)
    output:
        pca = "data/pca/all.poa.pcangsd.1dp4.cov"
    run:
        shell("gzip data/pca/all.poa.merged.1dp4.beagle")        
        shell("python tools/pcangsd/pcangsd.py -beagle data/pca/all.poa.merged.1dp4.beagle.gz -sites_save -o data/pca/all.poa.pcangsd.1dp4")


# (3b) run PCA with single read sampling approach
# -doIBS 1 prints a randomly sampled read from each individual at each position
# -doCov 1 prints out the covariance matrix which can be used for a PCA
# -minMaf 0.05 excludes sites with a minor allele freq less the 0.05
# -doMajorMinor 4 specifies the Major allele to be the reference allele specified by -ref
rule PCA_single:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/all_poa_bamlist.txt",
    output:
#        "data/pca/all.poa.1dp4.1e3.covMat"
        "data/pca/all.poa.1dp4.maf05.ibs.gz"
    params:
        prefix = "data/pca/all.poa.1dp4.maf05"
    run:
        shell("angsd \
        -GL 1 -P 15 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -setMinDepthInd 1 -setMaxDepthInd 4 \
        -minInd 8 \
        -bam {input.bamlist} \
        -doCounts 1 \
        -doMaf 1 \
        -minMaf 0.05 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}")
 #       -SNP_pval 1e-3 \
 #       -doCov 1 \
