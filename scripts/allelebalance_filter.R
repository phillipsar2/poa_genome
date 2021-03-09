library(dplyr)
library(tidyr)
library(splitstackshape)


# load table
## all poa
#snps <- read.csv('reports/filtering/all.poa.AB.table', header = T, na.strings=c("","NA"), sep = "\t") 
#samples = 8
## only P. pratensis
snps <- read.csv('reports/filtering/pPratensis.AB.table', header = T, na.strings=c("","NA"), sep = "\t") 
samples = 7

odd <- seq(3,2*samples+1,2)
even <- seq(4,2*samples+2,2)

# matrix of final read counts at each position
dp <- as.data.frame(matrix(nrow = dim(snps)[1], ncol = 3))
colnames(dp) <- c('CHROM', 'POS', 'AB')
dp$CHROM <- snps$CHROM
dp$POS <- snps$POS

# grab genotypes
gt <- snps[odd]

# identify hets
for (i in 1:dim(gt)[2]){
  gt[,i] <- gt[,i] %in% c('T/T', 'A/A', 'C/C', 'G/G')
}

# grab allele depth
ad <- snps[even]

# add factor level
for ( i in 1:dim(ad)[2]){
  levels(ad[,i]) <- c(levels(ad[,i]), "0,0")
}

# grab read counts for all invididuals
x <- colnames(ad)
all_read_counts <- cSplit(ad, x ,sep=',') %>% data.matrix()
dp$sum_format_dp <- rowSums(all_read_counts) # sum total filtered depth at each sites

# substitute zero reads for homozygotes
for (i in 1:dim(gt)[2]){
  for (l in 1:9999){
    if (gt[l,i] == TRUE){
      ad[l,i] <- '0,0'
    } 
  }
}

# Split ad columns by indiv
x <- colnames(ad)
read_counts <- cSplit(ad, x ,sep=',') # '_1' is the ref allele, '_2' is the alt allele

# Calculate ref and alt read counts
odd_reads <- seq(1,2*samples,2)
even_reads <- seq(2,2*samples+1,2)
  
ref_count <- read_counts[, ..odd_reads] %>% rowSums()
#alt_count <- rowSums(read_counts[,c(2,4,6,8,10,12,14,16)])
alt_count <- read_counts[, ..even_reads] %>% rowSums()

# Calculate AB
dp$AB <- alt_count/(ref_count + alt_count)

# Sum total reads from heterozygous indv
dp$total_het_read_count <- ref_count + alt_count




# ## Calculate total filtered read depth (sum of FORMAT/DP across all indiv)
# ad <- snps[even] # grab AD, again, keeping all individuals
# x <- colnames(ad)
# all_read_counts <- cSplit(ad, x ,sep=',')
# 
# odd_reads <- seq(1,2*samples,2) # sum reads across each row
# even_reads <- seq(2,2*samples+1,2)
# ref_count <- read_counts[, ..odd_reads] %>% rowSums()
# alt_count <- read_counts[, ..even_reads] %>% rowSums()
# 
# dp$sum_format_dp <- ref_count + alt_count
# 
# dp$sum_format_dp <- rowSums(all_read_counts)




# output table
## all poa
#write.table(dp, "reports/filtering/all.poa.AB.estimate.txt", row.names = FALSE)
## P. pratensis
write.table(dp, "reports/filtering/pPratensis.AB.estimate.txt", row.names = FALSE)
