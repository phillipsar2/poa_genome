# scripts provided by Silas

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

# generate angsd saf file (based on Silas's rule named pop_beagle)
rule angsd_saf:
    input:
        # reference, ancestral state ref, bams
        ref = config.ref,
        anc = config.ref,
        bams = "data/interm/mark_dups/bamlist.txt"
    output:
        # generating saf in parallel(?)
        saf = "data/angsd_pi/{popl}--{chrom}.saf.gz"
    params:
        scratch = "/scratch/aphillip/angsd",
        final = "data/angsd_pi/",
        # prefix for file names being generated
        prefix = "data/angsd_pi/{popl}--{chrom}",
        chrom = "{chrom}"
    shell:
        """
        mkdir -p {params.scratch}
#        rm -f {params.prefix}.arg
#        rm -f {params.prefix}.mafs.gz
#        rm -f {params.prefix}.saf.gz
#        rm -f {params.prefix}.saf.idx
#        rm -f {params.prefix}.saf.pos.gz

        module load angsd
        # samtools method for GL, request 15 threads,
        # unique reads only, remove reads with flag above 255, only proper pairs, don't trim, minimum mapping qual 40,minimum quality score 30
        # -doSaf = generate saf file, 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
        # -doMaf = estimate allele frequency
        # -doMajorMinor = how to decide the major allele, 4: for major alle according to reference states
        # -minInd = # individuals required to have data, 4 (50%) must have data at that site
        # -minMaf = , including singeltons can confound popu strucutre analysis (See https://www.biorxiv.org/content/10.1101/188623v2.full.pdf)
        angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -setMinDepth 1 -setMaxDepth 35 -minInd 4 -skipTriallelic 1 \
        -minMaf 0.1 \ 
        -ref {input.ref}  -anc {input.anc} \
        -doSaf 1 \
        -doMaf 2 \
        -doMajorMinor 4 \
        -r {params.chrom} \
        -bam {input.bams} -out {params.prefix}

#         mv {params.prefix}.arg {params.final}
#         mv {params.prefix}.mafs.gz {params.final}
#         mv {params.prefix}.saf.gz {params.final}
#         mv {params.prefix}.saf.idx {params.final}
#         mv {params.prefix}.saf.pos.gz {params.final}
        """

# generate SFS from saf file
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
## ? maybe change -minInd ?
rule angsd_pca:
    input:
        bams = "data/interm/mark_dups/bamlist.txt",
        ref = config.ref
    output:
        pca = "data/angsd_pi/pca/Ppratensis.pcangsd.cov"
    params:
        prefix_geno = "data/angsd_pi/Ppratensis",
        prefix_pca = "data/angsd_pi/pca/Ppratensis.pcangsd"
    run:
        shell("module load angsd")
        shell("angsd -GL 1 \
        -out {params.prefix_geno} \
        -nThreads 10 \
        -doGlf 2 -doMajorMinor 4 -SNP_pval 1e-6 -doMaf 2 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 -minInd 7 -skipTriallelic 1 \
        -bam {input.bams} \
        -ref {input.ref}")
        shell("python tools/pcangsd/pcangsd.py -beagle {params.prefix_geno}.beagle.gz -o {params.prefix_pca}")
