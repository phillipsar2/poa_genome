# scripts provided by Silas

# generate angsd saf file (based on Silas's rule named pop_beagle)
rule angsd_saf:
    input:
        # reference, ancestral state ref, bams
        ref = config.ref,
        anc = "data/anc/{ref}_anc.fa",
        bams = "data/interm/mark_dups/bamlist.txt"
    output:
        # generating saf in parallel(?)
        saf = "data/angsd_pi/{pop}--{chrom}.saf.gz"
    params:
        scratch = my_scratch,
#        final = "data/angsd_pi/",
        # prefix for file names being generated
        prefix = "data/angsd_pi/{pop}--{chrom}",
        chrom = "{chrom}"
    shell:
        """
#        mkdir -p {params.scratch}
#        rm -f {params.prefix}.arg
#        rm -f {params.prefix}.mafs.gz
#        rm -f {params.prefix}.saf.gz
#        rm -f {params.prefix}.saf.idx
#        rm -f {params.prefix}.saf.pos.gz

        module load angsd
        # samtools method for GL, request 15 threads,
        angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -ref {input.ref}  -anc {input.anc} \
        # -doSaf = generate saf file, 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
        -doSaf 1 \
        # -doMaf = Frequency (fixed major unknown minor)
#        -doMaf 2 \
        # -doMajorMinor = how to decide the major allele, 4: for major alle according to reference states 
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
        saf = "data/angsd_pi/{pop}--{chrom}.saf.gz"
    output:
        sfs = "data/angsd_pi/{pop}--{chrom}.sfs"
    params:
        prefix = "data/angsd_pi/{pop}--{chrom}",
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
        sfs = "data/angsd_pi/{pop}--{chrom}.sfs"
    output:
        pi = "data/angsd_pi/{pop}--{chrom}.{window}BP_theta.thetasWindow.gz.pestPG"
    params:
        prefix_in = "data/angsd_pi/{pop}--{chrom}",
        prefix_out = "data/angsd_pi/{pop}--{chrom}.{window}BP_theta",
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
