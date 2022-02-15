# Plot population structure and nucleotide diversity results 
# Author: Alyssa Phillips
# Last edited: 9/14/21

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(gridExtra)
library(nationalparkcolors)
library(paletteer)
library(gtable)

# > PCA ----

# >> Load metadata ----
pop<-read.table("~/Downloads/pop.info")
pop[,3] <- as.factor(c("P. compressa", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis", "P. pratensis"))
pop[,4] <- as.factor(c("Sorg", "Sorg", "Sorg", "Andro", "Schiz", "Schiz", "Andro", "Schiz"))
pop[,5] <- as.factor(c("1281-4", "1281-6", "1281-9", "1283-5", "1282-10", "1282-11", "CAM1428", "CAM1407:"))

# >> Load pca from PCangsd ----
## angsd
C <- as.matrix(read.table("~/poa/PCA/all.poa.pcangsd.1dp4.cov"))

# >> Compute eigenvalues and corresponding eigenvectors of S ----
# these are the principal components of S
e <- eigen(C)

# >> Print proportion of total var explained by the components ----
for (s in e$values) {
  print(s / sum(e$values))
}

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

# >> Plot with ggplot ----
df <- as.data.frame(cbind(e$vectors, pop))
names(df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "Population", "genotype", "Species", "bleh", "blah")
head(df)
str(df)

g1 <- ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = Population, fill = Species, color = Species),
             size = 3,
             alpha = 0.75,
             position = position_jitter(w = 0.02, h = 0.02)) +
  scale_shape_manual(values = c(24,23,21)) +
  scale_fill_manual(values = c("white","black")) +
  scale_color_manual(values = c("black", "white"))+
#  theme_ipsum(axis_title_just = "cc",
#              axis_title_size = 12,
#              grid_col = "black") + 
  
 # scale_fill_paletteer_d("nationalparkcolors::CraterLake", 1) +
#  scale_fill_grey(start = 1, end = 0) +
  xlab("PC1 (40.5%)") +
  ylab("PC2 (12.7%)") +
  ylim(-0.75,1) +
  xlim(-0.5,1) +
  theme_bw(base_size = 12) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

 
g1


# > Nucleotide diversity ----

# >> Load raw data ----
boulder = read.csv("~/poa/pi/boulder.1dp4.thetas.idx.pestPG", sep = "\t")
head(boulder)

crosspops = read.csv("~/poa/pi/all.crosspops.2dp6.thetas.idx.pestPG", sep = "\t")
head(crosspops)

# >> Discard windows that have nSites < 10% of the total window size ----
window = 1000
boulder$percent <- boulder$nSites / window
crosspops$percent <- crosspops$nSites / window

subdf_boulder <- boulder %>% filter(percent >= 0.1)
subdf_crosspops <- crosspops %>% filter(percent >= 0.1)

# >> Calculate pi while controlling for the number of sites with data in each window ----
# per site theta = (window theta) / (number of sites)

subdf_boulder$per_site <- subdf_boulder$tP / subdf_boulder$nSites
head(subdf_boulder$per_site)

subdf_crosspops$per_site = subdf_crosspops$tP/subdf_crosspops$nSites

# >> Calculate global pi ---- 
## Sum the per_site values and divide by the number of windows
boulder_global_pi <- sum(subdf_boulder$per_site)/dim(subdf_boulder)[1]
boulder_sd <- sd(subdf_boulder$per_site)

paste("boulder:", boulder_global_pi, "SD = ",boulder_sd )


crosspops_global_pi <- sum(subdf_crosspops$per_site)/dim(subdf_crosspops)[1]
crosspops_sd <- sd(subdf_crosspops$per_site)

paste("crosspops:", crosspops_global_pi, "SD = ", crosspops_sd)

means = data.frame(group = c("Boulder", "All populations"), per_site = c(boulder_global_pi, crosspops_global_pi))
means

head(subdf_boulder)

# >> Plot the per_site pi value for each window ----
pi <- as.data.frame(subdf_boulder$per_site)
pi$group <- "Boulder"
colnames(pi) <- c("per_site", "group")
pi2 <- as.data.frame(subdf_crosspops$per_site)
pi2$group <- "All populations"
colnames(pi2) <- c("per_site", "group")

yum <- rbind(pi, pi2)
head(yum)

g2 <- ggplot(yum, 
             aes(x=group, 
                 y = per_site,
                 fill = "black"
             )) +
                 # fill = group,
                 # color = group)) +
  geom_violin(
    fill = "darkgray"
  ) +
  ylim(0, 0.04) +
  # scale_color_paletteer_d("nationalparkcolors::CraterLake") +
  # scale_fill_paletteer_d("nationalparkcolors::CraterLake") +
  # theme_ipsum(axis_title_just = "cc",
  #             axis_title_size = 12,
  #             grid_col = "black",
  #             axis_col = "black") +
  theme_bw(base_size = 12) +
  theme(legend.position="none",
        plot.margin = unit(c(10,10,0,10), "pt"),
        axis.title.y = element_blank()
        ) +
  ylab("Per-site nucleotide diversity") +
  coord_flip()+
  geom_point(data = means, col = 'black', size = 3)

g2

## > Combine plots into one figure ----
gg <- grid.arrange(g1, g2, ncol = 1,
             heights = c(2.5,1))

# > Save figure ----
ggsave("~/poa/Paper/pi_pca_fig_grayscale_11042021.jpg", plot = gg, dpi = 300, width = 5.25, height =5, units = "in")

