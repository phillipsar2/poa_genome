import config

# Sample names
SAMPLE = glob_wildcards("data/raw/sequences/{sample}_1.fq.gz").sample
#print(SAMPLE)

# Treating poa reference genome as the sample
#SAMPLE = ["poa-v1"]
#print(SAMPLE)


# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

# Set chromosome IDs
CHROM = ["Super-Scaffold_107", "Super-Scaffold_475", "Super-Scaffold_13795", "Super-Scaffold_40012", "Super-Scaffold_1000001", "Super-Scaffold_1000002", "Super-Scaffold_1000003", "Super-Scaffold_1000004", "Super-Scaffold_1000005", "Super-Scaffold_1000006", "Super-Scaffold_1000007", "Super-Scaffold_1000008", "Super-Scaffold_1000009", "Super-Scaffold_1000011", "Super-Scaffold_1000012", "Super-Scaffold_1000013", "Super-Scaffold_1000014", "Super-Scaffold_1000015", "Super-Scaffold_1000016"]

# Set window size for angsd sliding window statistics
#WINDOW = ["50000"]

# Set population. Options: ["all", "boulder", "manitoba"]
POPL = ["all"]

rule all:
    input:
        # Aligning reads
#        markdups = expand(config.mark_dups, sample = SAMPLE),
        # Assess quality of mapped reads
#        bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLE),
        # Call SNPs
#         mpileup = "data/bcf/all.poa.raw.bcf",
#         bcf2vcf = "data/vcf/mpileup/all.poa.raw.vcf",
#         diag = "reports/filtering/all.poa.snps.vcf.table",       
#        haplo_caller = expand("data/vcf/{sample}.vcf", sample = SAMPLE),
#        split_int = expand("data/processed/scattered_intervals/{interval}-scattered.interval_list", interval = INTERVALS),
#        joint_geno = expand("data/raw/vcf_bpres/{interval}.raw.vcf", interval = INTERVALS),
        # Filter SNPs
#         hard_filt = "data/processed/filtered_snps/all.poa.filtered.nocall.snps.vcf",
#         diag_depth = "reports/filtering/all.poa.depth.filtered.nocall.table",
#         filter_depth = "data/processed/filtered_snps/all.poa.filtered.nocall.2dp20.snps.vcf",
#        dp_nocall = "data/processed/filtered_snps/all.poa.filtered.nocall.2dp20.max0.snps.vcf",
#        grab_Ppr = "data/processed/filtered_snps/poa.pratensis.filtered.nocall.2dp20.max0.snps.vcf",
#        diag_AB = config.ab_table,
#        calc_AB = "reports/filtering/all.poa.AB.estimate.txt",
        # Analysis of pop panel with angsd
#        to_beagle = expand(config.beagle, chrom = CHROM),
#        pca = "data/angsd_pi/pca/Ppratensis.pcangsd.cov"
#        sfs = expand(config.sfs, chrom = CHROM),
#        pi = expand(config.stats, chrom = CHROM),
##### Rules for aligning pacbio reads
#        mpileup = expand("data/vcf/pacbio.{chrom}.vcf", chrom = CHROM),
#        diag = "reports/bamqc/pacbio/qualimapReport.html",
##### Rules for identifying the APOSTART gene
#        bwa = "data/interm/mapped_bam/apostat.mapped.bam",
##### Rules for identifying genes and building consensus sequences ####
       # Call vcf
        haplo = expand("data/vcf/trnLtrnF/{sample}.trnLtrnF.vcf", sample = SAMPLE),
        # SNP calling
#        hap_vcf = expand("data/vcf/ETS/{sample}.vcf", sample = SAMPLE),
#        hard_filt = expand("data/processed/filtered_snps_bpres/{sample}.ETS.filtered.snps.vcf", sample = SAMPLE),
        # Build consensus seq
#        alt_ref = expand("data/processed/pseudo_ref/{sample}.pseudo.fasta", sample = SAMPLE)
#        ETS = expand("data/gene/ETS/{sample}.consensus.ETS.fasta",sample = SAMPLE)
        # Align gene to consensus seq
#        bwa_single = expand("data/interm/mapped_bam/{sample}.rpoB-trn-C.mapped.bam", sample = SAMPLE)
#        get_seq = expand("data/processed/gene/trnTtrnLtrnF/{sample}.trnTtrnLtrnF.fasta", sample = SAMPLE)
        # Multiple sequence alignment
#        maaft = config.multi_aln

        # Prepping vcf for consensus sequence
#        vcf = expand("data/processed/filt_indels/{sample}.vcf.gz", sample = SAMPLE),
#        txt = expand("data/processed/iupac/{sample}.txt", sample = SAMPLE),
#        head = expand("data/processed/filt_indels/{sample}.header", sample = SAMPLE)

#        dp = expand("data/processed/filtered_snps_bpres/{intervals}.filtered.dp3_77.snps.vcf", intervals = INTERVALS),
#        dp2 = expand("data/processed/filtered_snps_bpres/{intervals}.filtered.dp3_77.nocall.snps.vcf", intervals = INTERVALS),
#        bgzip_vcf = expand("data/processed/filtered_snps_bpres/{intervals}.filtered.dp3_77.nocall.snps.vcf.gz", intervals = INTERVALS),
#        combine = "data/processed/filtered_snps_bpres/oryza_glum.vcf.gz",
#        depth_diag = "reports/filtering_bpres/oryza_glum.table"
#        wholegenome = "data/processed/filtered_snps_bpres/oglum_wholegenome.vcf.gz"

# Rules
#include: "rules/mapping.smk"
include: "rules/process_bam.smk"
#include: "rules/calling.smk"
include: "rules/consensus_seq.smk"
#include: "rules/filtering.smk"
#include: "rules/multiseq_align.smk"
#include: "rules/pop_gen.smk"
#include: "rules/ref_AB.smk"
#include:"rules/apostart.smk"
