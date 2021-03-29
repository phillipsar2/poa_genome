### Calling SNPs for whole genome ###

# Genome
#ref = "data/genome/poa_01062021/Ppratensis.normalized.nonull.fasta"
#contig_list = "data/genome/poa_01062021/Ppratensis.contig.list"
#sample_map = "data/processed/poa.sample_map"

# Process bams
## Process population panel
#sort_in = "data/interm/mapped_bam/{sample}.mapped.bam"
#sort_out = "data/interm/sorted_bam/{sample}.sorted.bam"
#add_rg = "data/interm/addrg/{sample}.rg.bam"
#mark_dups = "data/interm/mark_dups/{sample}.dedup.bam"

### delete this section
## Process the pacbio reads
#sort_in = "data/interm/mapped_bam/pacbio.mapped.bam"
#sort_out = "data/interm/sorted_bam/pacbio.sorted.bam"
#add_rg = "data/interm/addrg/pacbio.rg.bam"
#mark_dups = "data/interm/mark_dups/pacbio.dedup.bam"
### delete this section


# Evaluating allele balance and depth - all Poa
#ab_table = "reports/filtering/all.poa.AB.table"

# Evaluating allele balance and depth - Ppratensis
ab_table = "reports/filtering/pPratensis.AB.table"


# Pop gen - Ppratensis
## convert vcf to beagle
final_vcf = "data/processed/filtered_snps/poa.pratensis.filtered.nocall.2dp20.max0.snps.vcf"
beagle = "data/beagle/ppratensis/poa.pratensis.{chrom}.beagle.gz"

# Pop gen - all Poa
## convert vcf to beagle
#final_vcf = "data/processed/filtered_snps/all.poa.filtered.nocall.2dp20.max0.snps.vcf"
#beagle = "data/beagle/all-poa/all.poa.{chrom}.BEAGLE.PL"

### Consensus sequence building ###

# Pulling out the chloroplast pseudogenomes
#sort_out = "data/sorted_bam/{sample}.TLF.ref.sorted.bam"

# bwa_pull - aligning individual's reads to the reference sequence
gene = "data/gene/trnTtrnLtrnF/DQ354006.1.fasta" # reference
sort_in = "data/interm/mapped_bam/{sample}.trnLtrnF.mapped.bam"

sort_out = "data/interm/sorted_bam/{sample}.trnLtrnF.sorted.bam"
add_rg = "data/interm/addrg/{sample}.trnLtrnF.rg.bam"
mark_dups = "data/interm/mark_dups/{sample}.trnLtrnF.dedup.bam"

haplo = "data/vcf/trnLtrnF/{sample}.trnLtrnF.vcf"

######### Calculate theta ###############

# Boulder population
#bamlist = "data/interm/mark_dups/boulder_bamlist.txt"
#saf = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.saf.gz"
#prefix = "data/angsd_pi/boulder/pratensis.boulder.{chrom}"
#ind = 5
#sfs = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.sfs"
#stats = "data/angsd_pi/boulder/pratensis.boulder.{chrom}.thetas.idx.pestPG"

# Across populations
#bamlist = "data/interm/mark_dups/crosspops_bamlist.txt"
#saf = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.saf.gz"
#prefix = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}"
#ind = 3
#sfs = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.sfs"
#stats = "data/angsd_pi/crosspops/pratensis.crosspops.{chrom}.thetas.idx.pestPG"
