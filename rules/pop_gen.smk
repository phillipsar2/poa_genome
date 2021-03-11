# scripts modified from Silas Tittes's shared scripts

# convert filtered vcf to beagle format
rule vcf_to_beagle:
    input:
        vcf = config.final_vcf
    output:
        beagle = config.beagle
    params:
        chr = "{chrom}"
    run:
        shell("module load vcftools")
        shell("vcftools --vcf {input.vcf} --BEAGLE-PL --stdout --chr {params.chr} > {output}")

# subset Boulder individuals
rule grab_boulder:
    input:
        vcf = "data/proccessed/filtered_snps/poa.pratensis.filtered.nocall.2dp20.max0.snps.vcf",
        samples = "data/raw/boulder_samples.txt"
    output:
        "data/beagle/boulder/pratensis.boulder.{CHROM}.BEAGLE.PL.gz"
    params:
        chr = "{chrom}"
    run:
        shell("bcftools view -S {input.samples} {input.vcf} | \
        vcftools --vcf - --BEAGLE-PL --stdout --chr {params.chr} > {output}")

# subset Manitoba individuals + 1 Boulder individual
rule grab_crosspops:
    input:
        vcf = "poa.pratensis.filtered.nocall.2dp20.max0.snps.vcf",
        samples = "data/raw/crosspops_samples.txt"
    output:
        "data/beagle/crosspops/pratensis.crosspops.{CHROM}.BEAGLE.PL.gz"
    params:
        chr = "{chrom}"
    run:
        shell("bcftools view -S {input.samples} {input.vcf} | \
        vcftools --vcf - --BEAGLE-PL --stdout --chr {params.chr} > {output}")


## Generating saf file - site allele frequency likelihood
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
#        beagle = "data/beagle/boulder/pratensis.boulder.merged.BEAGLE.PL.gz"
        ref = config.ref,
        bamlist = config.bamlist
    output:
        config.saf
    params:
        ind = {config.ind},
        prefix = config.prefix,
        chrom = "{chrom}"
    run:
#        shell("angsd -doSaf 4 \
#        -beagle {input.beagle} \     
#        -out {params.prefix}")
        shell("angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 \
        -setMinDepthInd 2 -setMaxDepthInd 20 \ 
        -minInd {params.ind} \  
        -ref {input.ref} -anc {input.ref} \ 
        -doSaf 1 \
        -r {params.chrom}
        -bam {input.bamlist} \ 
        -out {params.prefix}")


## Generate SFS (site frequency spectrum) from saf file
rule pop_sfs:
    input:
        ref = config.ref,
        bams = "data/interm/mark_dups/bamlist.txt",
        saf = "data/angsd_pi/{popl}--{chrom}.saf.gz"
    output:
        sfs = "data/angsd_pi/{popl}--{chrom}.sfs"
    params:
        prefix = "data/angsd_pi/{popl}--{chrom}",
    shell:
        """
        module load angsd
        # generate sfs from saf, 10 threads, maximum likelihood estimate, unfolded
        realSFS {params.prefix}.saf.idx  -P 10 > {params.prefix}.sfs
        # calculate thetas for each site
        realSFS saf2theta {params.prefix}.saf.idx -sfs {params.prefix}.sfs -outname {params.prefix}
        """

# calculate pi
rule pop_pi:
    input:
        sfs = "data/angsd_pi/{popl}--{chrom}.sfs"
    output:
        pi = "data/angsd_pi/{popl}--{chrom}.{window}BP_theta.thetasWindow.gz.pestPG"
    params:
        prefix_in = "data/angsd_pi/{popl}--{chrom}",
        prefix_out = "data/angsd_pi/{popl}--{chrom}.{window}BP_theta",
        win = "{window}"
    shell:
        """
        module load angsd
        # calculate Tajimas D and other stats in sliding windows
        thetaStat do_stat {params.prefix_in}.thetas.idx \
        -win {params.win} \
        -step {params.win} \
        -outnames {params.prefix_out}.thetasWindow.gz
        """

# run PCA
rule angsd_pca:
    input:
        beagle ="data/angsd_pi/Ppratensis.beagle.gz"
    output:
        pca = "data/angsd_pi/pca/Ppratensis.pcangsd.cov"
    run:
        shell("python tools/pcangsd/pcangsd.py -beagle {input.beagle} -o {output._pca}")
