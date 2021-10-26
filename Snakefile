import config

# Sample names
SAMPLE = glob_wildcards("data/raw/sequences/{sample}_1.fq.gz").sample
#print(SAMPLE)

# Gene names - specify one at a time
# options = ['ETS', 'ITS', 'matk', 'rpoB-trn-C', 'trnTtrnLtrnF']
GENE = ETS

# Treating poa reference genome as the sample
#SAMPLE = ["poa-v1"]
#print(SAMPLE)


# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

# Set chromosome IDs
CHROM = ["Super-Scaffold_107", "Super-Scaffold_475", "Super-Scaffold_13795", "Super-Scaffold_40012", "Super-Scaffold_1000001", "Super-Scaffold_1000002", "Super-Scaffold_1000003", "Super-Scaffold_1000004", "Super-Scaffold_1000005", "Super-Scaffold_1000006", "Super-Scaffold_1000007", "Super-Scaffold_1000008", "Super-Scaffold_1000009", "Super-Scaffold_1000011", "Super-Scaffold_1000012", "Super-Scaffold_1000013", "Super-Scaffold_1000014", "Super-Scaffold_1000015", "Super-Scaffold_1000016"]

# Set population. Options: ["all", "boulder", "manitoba"]
POPL = ["all"]

rule all:
    input:
### Genetic diversity rules
        # Aligning reads
        markdups = expand(config.mark_dups, sample = SAMPLE),
        # Assess quality of mapped reads
        bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLE),
        # Nucleotide diversity
#        sfs = expand(config.sfs, chrom = CHROM),
#        pi = expand(config.stats, chrom = CHROM),
        # PCA
#        beagle = expand("data/pca/all.poa.{chrom}.beagle.gz", chrom = CHROM),
#        pca = "data/pca/all.poa.pcangsd.2dp6.cov"
### Rules for calculating AB in the reference genome
#        mpileup = expand("data/vcf/pacbio.{chrom}.vcf", chrom = CHROM),
#        diag = "reports/bamqc/pacbio/qualimapReport.html",
#         merge_vcf = "pacbio.all_scaff.vcf.gz",
#        diag_ab = "reports/filtering/pacbio.AB.table",
#        calc_ab = "reports/filtering/pacbio.AB.estimate.txt",
### Rules for generating consensus marker sequences
       # Align to gene reference
#        map = expand(sort_in, sample = SAMPLE, gene = GENE)
       # Call SNPs
#        haplo = expand("data/vcf/{gene}/{sample}.{gene}.snps.vcf", sample = SAMPLE,  gene = GENE),
       # Filter snps
#        hard_filt = expand("data/processed/filtered_snps/{gene}/{sample}.{gene}.filtered.nocall.snps.vcf", sample = SAMPLE, gene = GENE)
        # Call consensus seq
#        consensus = expand("data/gene/{gene}/{sample}.consensus.{gene}.fasta", sample = SAMPLE, gene = GENE)

# Rules
include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/calling.smk"
#include: "rules/consensus_seq.smk"
#include: "rules/filtering.smk"
#include: "rules/pop_gen.smk"
#include: "rules/ref_AB.smk"
