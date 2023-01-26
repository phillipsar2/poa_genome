### Title: Population structure and nucleotide diversity results 
### Author: Alyssa Phillips
### Last edited: 9/14/21

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(gridExtra)
# library(nationalparkcolors)
library(paletteer)
library(gtable)
library(matrixStats)

###
### PCA ----
###

# > Load metadata ----
pop <- read.table("~/Downloads/pop.info")
pop[,3] <- as.factor(c("P. compressa", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis"))
pop[,4] <- as.factor(c("Sorg", "Sorg", "Sorg", "Andro", "Schiz", "Schiz", "Andro", "Schiz"))
pop[,5] <- as.factor(c("1281-4", "1281-6", "1281-9", "1283-5", "1282-10", "1282-11", "CAM1428", "CAM1407"))

# >> Subset to pratensis only ----
pop <- pop[pop$V3 == "P. pratensis",]
pop

# > Load ANGSD covariance matrix ----
# C <- as.matrix(read.table("~/poa/PCA/all.poa.1dp4.1e4.covMat")) # all poa
C <- as.matrix(read.table("~/poa/PCA/pratensis.1dp4.maf05.covMat")) # pratensis only

# > Compute eigenvalues and corresponding eigenvectors of S ----
# these are the principal components of S
e <- eigen(C)

# > Print proportion of total var explained by the components ----
pc <- e$values / sum(e$values)
pc

#pdf("all--PCAngsd1.pdf")
# 
# plot(e$vectors[,1:2],
#      col="black",
#      pch=c(23, 21)[as.numeric(pop[,3])],
#      xlab="PC1 (29.2%)",
#      ylab="PC2 (15.8%)",
#      cex = 2,
#      bg=c(0,3,4)[as.numeric(pop[,1])])
# #     main = "All poa, Super-Scaffold_1000016 ")
# legend("top", title="Location         Species", 
#        legend=c("Boulder", "Tolstoi, Manitoba", "Argyle, Manitoba","P. pratensis", "P. compressa"),
#        col=c(3,4,"black",rep("black",2)),
#        pch=c(21,21,21,23, 21),
#        bty="n", border=F, ncol=2)

# > Plot with ggplot ----
df <- as.data.frame(cbind(e$vectors, pop))
# names(df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "Population", "genotype", "Species", "bleh", "blah")
names(df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Population", "genotype", "Species", "bleh", "blah")
head(df)
str(df)

g1 <- ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Population, 
                 # fill = Species,
                 # color = Species
                 ),
             size = 4,
             alpha = 0.75,
             ) +
  # scale_shape_manual(values = c(24,23,21)) + # all poa
  scale_shape_manual(values = c(17,18,19)) + # pratensis only
  # scale_fill_manual(values = c("white","black")) +
  # scale_color_manual(values = c("black", "white"))+
#  theme_ipsum(axis_title_just = "cc",
#              axis_title_size = 12,
#              grid_col = "black") + 
  
 # scale_fill_paletteer_d("nationalparkcolors::CraterLake", 1) +
#  scale_fill_grey(start = 1, end = 0) +
  xlab(paste0("PC1 (", round(pc[1]*100,1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2]*100, 1), "%)")) +
  ylim(-0.75,1) +
  xlim(-0.5,1) +
  # labs(tag = "A") +
  theme_bw(base_size = 12) +
  theme(
    legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    plot.tag = element_text())
 
g1

ggsave("~/poa/Paper/prat_pca_01232023.jpeg")

###
###  Nucleotide diversity per genome ----
###

# > Load the genotype table ----
# There should be no missing data (-1 value). 
r1 <- read.table(gzfile("~/poa/pi/pratensis.1dp4.MQ20.r1.ibs.gz"), header= T) 
r2 <- read.table(gzfile("~/poa/pi/pratensis.1dp4.MQ20.r2.ibs.gz"), header= T) 
r3 <- read.table(gzfile("~/poa/pi/pratensis.1dp4.MQ20.r3.ibs.gz"), header= T) 

dim(r1) # number of SNPs

# > Calculate pi ----
# Sum columns (each individual) to get pi per genome and divide by number of snps
# Remove rows with missing data
# pi_ind <- colSums(g_nomiss[,5:dim(g_nomiss)[2]]) / dim(g_nomiss)[1]

calc_pi <- function(g){
  g[g == -1 ] <- NA
  g_nomiss <- g[complete.cases(g),]
  pi <- colSums(g_nomiss[,5:dim(g_nomiss)[2]]) / dim(g_nomiss)[1]
  return(pi)
}

r1_pi <- calc_pi(r1)
r2_pi <- calc_pi(r2)
r3_pi <- calc_pi(r3)

pi_ind_df <- cbind(pop[2:8,], r1_pi, r2_pi, r3_pi)
colnames(pi_ind_df) <- c("pop", "geno", "species", "contaminant", "ID", "pi_run1", "pi_run2", "pi_run3")
pi_ind_df$metapop <- c("Boulder","Boulder","Boulder","Boulder","Boulder","Manitoba","Manitoba")
pi_ind_df

# Test if the two pops are different?
# pratensis <- pi_df %>%
  # filter(species == "P. pratensis" )
# t.test(pi ~ metapop, data = pi_ind_df)
'The 95% confidence interval is very large and crosses zero'

# > Variance across runs ----
rowVars(as.matrix(pi_ind_df[,6:8])) %>% mean()

# > Summary stats ----
# Overall mean
rowMeans(as.matrix(pi_ind_df[,6:8])) %>% mean()

# Overall sd
stack(pi_ind_df[,6:8])[,1] %>% sd()

# Within population ranges
pi_ind_df[pi_ind_df$pop == "Boulder", 6:8] %>% 
  as.matrix() %>% 
  rowMeans() %>%
  range()

pi_ind_df[pi_ind_df$metapop == "Manitoba", 6:8] %>% 
  as.matrix() %>% 
  rowMeans() %>%
  range()

# > Plot ----
g2 <- ggplot(pi_ind_df) +
  geom_segment( aes(x = geno, xend = geno, y = pi_run1, yend = pi_run2), color = "black") +
  geom_segment( aes(x = geno, xend = geno, y = pi_run1, yend = pi_run3), color = "black") +
  geom_point( aes(x = geno, y = pi_run1), color = "black", size = 1 ) +
  geom_point( aes(x = geno, y = pi_run2), color = "black", size = 1 ) +
  geom_point( aes(x = geno, y = pi_run3), color = "black", size = 1 ) +
  facet_wrap(~metapop, ncol = 1, scale = "free") +
  coord_flip()+
  theme_bw(base_size = 12) +
  labs(tag = "B") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(10,10,0,10), "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text()
  ) +
  ylim(c(0.0064, 0.007)) +
  ylab("Nucleotide diversity per genome")

g2

# pi_ind_df %>%
#   # filter(species == "P. pratensis" ) %>%
#   ggplot(aes(x = metapop, y = pi)) +
#     # geom_boxplot() +
#     # geom_jitter(width = 0.1) +
#     geom_point() +
#     theme_bw(base_size = 12) +
#     theme(legend.position="none",
#         plot.margin = unit(c(10,10,0,10), "pt"),
#         axis.title.x = element_blank()) +
#   ylab ("Nucleotide diversity per genome")

# > Plot the per_site pi value for each window 
# pi <- as.data.frame(subdf_boulder$per_site)
# pi$group <- "Boulder"
# colnames(pi) <- c("per_site", "group")
# pi2 <- as.data.frame(subdf_crosspops$per_site)
# pi2$group <- "All populations"
# colnames(pi2) <- c("per_site", "group")
# 
# yum <- rbind(pi, pi2)
# head(yum)
# 
# g2 <- ggplot(yum, 
#              aes(x=group, 
#                  y = per_site,
#                  fill = "black"
#              )) +
#                  # fill = group,
#                  # color = group)) +
#   geom_violin(
#     fill = "darkgray"
#   ) +
#   ylim(0, 0.04) +
#   theme_bw(base_size = 12) +
#   theme(legend.position="none",
#         plot.margin = unit(c(10,10,0,10), "pt"),
#         axis.title.y = element_blank()
#         ) +
#   ylab("Per-site nucleotide diversity") +
#   coord_flip()+
#   geom_point(data = means, col = 'black', size = 3)
# 
# g2

###
### Combine plots into one figure ----
###
gg <- grid.arrange(g1, g2, ncol = 1,
             heights = c(4,3))

# > Save figure ----
ggsave("~/poa/Paper/pi_pca_fig_grayscale_01182023.jpg", plot = gg, dpi = 300, width = 5.5, height =8, units = "in")
