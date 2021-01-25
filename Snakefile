import config

# Sample names
SAMPLE = glob_wildcards("data/raw/sequences/{sample}_1.fq.gz").sample
#print(SAMPLE)

# Treating poa reference genome as the sample
#SAMPLE = ["poa-v1"]
#print(SAMPLE)

# Set number of intervals for gatk to 200
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

rule all:
    input:
        # Aligning reads
        markdups = expand(config.mark_dups, sample = SAMPLE),
        # Assess quality of mapped reads
#        bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLE),
        # Call SNPS
        haplo_caller = expand("data/vcf/{sample}.vcf", sample = SAMPLE),
#        split_int = expand("data/processed/scattered_intervals/{interval}-scattered.interval_list", interval = INTERVALS),
##### Rules for identifying genes and building consensus sequences ####
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
include: "rules/mapping.smk"
include: "rules/process_bam.smk"
include: "rules/calling.smk"
#include: "rules/consensus_seq.smk"
#include: "rules/filtering.smk"
#include: "rules/multiseq_align.smk"
