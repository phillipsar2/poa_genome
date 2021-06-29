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


## Nucleotide diversity ----

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
        -setMinDepthInd 2 -setMaxDepthInd 6 \
        -minInd {params.ind} \
        -ref {input.ref} -anc {input.ref} \
        -doSaf 1 \
        -r {params.chrom} \
        -bam {input.bamlist} \
        -out {params.prefix}")


## Generate global estimate of SFS (site frequency spectrum) from saf file
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

# calculate thetas (and neutrality tests) in 10k sliding windows
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

### PCA --------

# generate GL in beagle format
# doMajorMinor 4: use refence allele as major
# assuming a mean coverage of 2X, 99% of sites should have 6X coverage in a poission distribution; ppois(6, lambda = 2)
rule angsd_beagle:
    input:
        ref = config.ref,
        bamlist = "data/interm/mark_dups/bamlist.txt"
    output:
        "data/pca/all.poa.{chrom}.2dp6.beagle.gz"
    params:
        prefix = "data/pca/all.poa.{chrom}.2dp6",
        chrom = "{chrom}"
    run:
        shell("angsd -GL 1 -P 15 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 \
        -setMinDepthInd 2 -setMaxDepthInd 6 \
        -minInd 8 \
        -ref {input.ref} -anc {input.ref} \
        -r {params.chrom} \
        -bam {input.bamlist} \
        -out {params.prefix}")

rule angsd_pca:
    input:
        beagle = expand("data/pca/all.poa.{chrom}.2dp6.beagle.gz", chrom = CHROM)
    output:
        pca = "data/pca/all.poa.pcangsd.2dp6.cov"
    run:
#        shell("cat <(zcat data/pca/*beagle.gz | head -n1) <(zcat data/pca/*beagle.gz | \
#        grep -v marker) > data/pca/all.poa.merged.2dp6.beagle")
        shell("gzip data/pca/all.poa.merged.2dp6.beagle")        
        shell("python tools/pcangsd/pcangsd.py -beagle data/pca/all.poa.merged.2dp6.beagle.gz -o data/pca/all.poa.pcangsd.2dp6")
