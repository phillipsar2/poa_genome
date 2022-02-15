# Building tree figures from tree files for poa genome project
# author: Alyssa Phillips
# date: May 14, 2021
# requires 32G memory
# R/4.0.2

library(ape)
#BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)

##https://yulab-smu.top/treedata-book/chapter4.html#displaying-tree-components

# > Load tree with treeio ----
# NEXUS format from Mr.Bayes
# output parsing supported by FigTree and treeio
bayes_tree <- read.mrbayes("~/poa/Trees/ITS/Poa_genome_ITS_6Apr2021.tre")

bayes_tree

# return vector of features/attirbutes
get.fields(bayes_tree)

# get features/attributes
get.data(bayes_tree)

# > Edit tree file ----
# >> Convert into tidy object ----
x <- as_tibble(bayes_tree)

# find rows that need their names changed
grepl("v1", x$label)
grepl("AN18", x$label)

# >> Alter sequence names  ----

# ITS tree
x[32,4]<- "Reference genome"
x[1,4] <- "AN18N065"
x[25,4] <- "AN18N067"
x[26,4] <- "AN18N070"
x[27,4] <- "AN18N123"
x[28,4] <- "AN18N353"
x[29,4] <- "AN18N107"
x[30,4] = "AN18N124"
x[31,4] = "AN18N234"

# # TLF tree
# x[48,4]<- "Reference genome"
# x[70,4] <- "AN18N065"
# x[53,4] <- "AN18N067"
# x[54,4] <- "AN18N070"
# x[55,4] <- "AN18N123"
# x[59,4] <- "AN18N353"
# x[58,4] <- "AN18N107"
# x[56,4] = "AN18N124"
# x[57,4] = "AN18N234"

# # ETS tree
# x[18,4]<- "Reference genome"
# x[89,4] <- "AN18N065"
# x[8,4] <- "AN18N067"
# x[9,4] <- "AN18N070"
# x[11,4] <- "AN18N123"
# x[13,4] <- "AN18N353"
# x[10,4] <- "AN18N107"
# x[12,4] = "AN18N124"
# x[14,4] = "AN18N234"



# >> Convert back to treedata (not phylo) object ----
mytree <- as.treedata(x)
#mytree <- as.phylo(x) # will lose branch support values

# >> Drop AN18N112 tip ----
# low coverage, not included in final data set

#final_tree <- treeio::drop.tip(mytree, "AN18N112_ITS")
#final_tree <- treeio::drop.tip(mytree, "AN18N112_trnTtrnLtrnF")

# >> Generate metadata ----

# TLF metadata
tip_metadata <- matrix(nrow = length(x$label), ncol = 2)
tip_metadata[,1] <- x$label
tip_metadata[32,2] = "Y"
tip_metadata[1,2]<- "Y"
tip_metadata[25,2] <- "Y"
tip_metadata[26,2] <- "Y"
tip_metadata[27,2] <- "Y"
tip_metadata[28,2] <- "Y"
tip_metadata[29,2] <- "Y"
tip_metadata[30,2] = "Y"
tip_metadata[31,2] = "Y"

# # ETS metadata
# tip_metadata <- matrix(nrow = length(x$label), ncol = 2)
# tip_metadata[,1] <- x$label
# tip_metadata[18,2] = "Y"
# tip_metadata[89,2]<- "Y"
# tip_metadata[8,2] <- "Y"
# tip_metadata[9,2] <- "Y"
# tip_metadata[11,2] <- "Y"
# tip_metadata[13,2] <- "Y"
# tip_metadata[10,2] <- "Y"
# tip_metadata[12,2] = "Y"
# tip_metadata[14,2] = "Y"

# # ITS metadata
# tip_metadata <- matrix(nrow = length(x$label), ncol = 2)
# tip_metadata[,1] <- x$label
# tip_metadata[2,2] = "Y"
# tip_metadata[62,2]<- "Y"
# tip_metadata[47,2] <- "Y"
# tip_metadata[48,2] <- "Y"
# tip_metadata[49,2] <- "Y"
# tip_metadata[50,2] <- "Y"
# tip_metadata[51,2] <- "Y"
# tip_metadata[52,2] = "Y"
# tip_metadata[53,2] = "Y"

tip_metadata = as.data.frame(tip_metadata)
colnames(tip_metadata) <- c("taxa", "poa_sample")

# > Plot full tree ---- 

# Reroot tree to Poa_wolfii_G1065_S5800 (node 160)
final_tree <- reroot(mytree, 160)

t <- ggtree(final_tree) +
  ## ETS + ITS: x=0.18, y=0
  ## TLF: x=0.04, y=0
  geom_treescale(x=0.18, y=0, fontsize = 2) + # scale branch lengths
  geom_tiplab(size = 2) +
  ## ITS: 0, 0.25
  ## ETS: 0.2
  ## TLF: 0.5
  ggplot2::xlim(0, 0.25) +
  # geom_nodelab(aes(x=branch, label=round(as.numeric(prob),2)),
  #             size = 1.5, vjust = -.5)
  geom_nodelab(aes(x=branch, label=node),
             size = 1.5, vjust = -.5)

# Add meta data to tree
t %<+% tip_metadata + 
  geom_tippoint(aes(color=poa_sample), size=1) + 
  theme(legend.position="NA") +
  scale_color_manual(values = c("deepskyblue", "black")) +
  # Highlight clades
  ## ITS: 153 & 157, extend = 0.08, extend = 0.09
  ## ETS: 139 & 176, extend = 0.08, extend = 0.09
  ## TLF: 152 & 166, extend = 0.012, extend = 0.01
  geom_hilight(node = 141, fill = "darkgray", alpha = 0.3, extend = 0.08) +
  geom_hilight(node = 130, fill = "darkgray", alpha = 0.3, extend = 0.06)

ggsave("~/poa/Trees/ITS.tree.11012021.jpeg", width = 8.5, height = 11, dpi = 300, units = "in" )
    




# > Small tree for figure ----

# >> List of species to keep in small tree ----

# select branches to keep in subsetted tree
# final_tree$tip.label %>% sort()

# species <- c("AN18N067","AN18N070","AN18N123","AN18N353","AN18N107","AN18N124", "AN18N234","Reference genome", 
#              "Poa_pratensis_angustifolia_G0797_Catalan03_16",
#              "Poa_pratensis_G1246_S7499",
#              "Poa_pratensis_pratensis_G0976_Olonova0342",
#              "Poa_pratensis_pratensis_G0170_6291",
#              "Poa_pratensis_pratensis_G0178_6310",
#              "Poa_pratensis_pratensis_G2142_10592",
#              "Poa_pratensis_irrigata_G1119_S6044",
#              "Poa_planifolia_G1011_Peterson19233",
#              "Poa_costiniana_G0615_7356_1",
#              "Poa_poiformis_G0623_7381",
#              "Poa_abbreviata_G0046_5816",
#              "Poa_porsildii_G1157_S6147_1" ,
#              "Poa_sibirica_sibirica_G0977_Olonova2003_45" ,
#              "Poa_remota_G1254_S7540",
#              "Poa_arachnifera_G1066_S5801" ,
#              "Poa_fendleriana_G0171_6292",
#              "Poa_occidentalis_G1382_Peterson18918" ,
#              "Poa_pratensis_alpigena_G0036_5801",
#              "Poa_reflexa_G1229_S7422",
#              "Poa_chaixii_G1049_S4677",
#              "Poa_pseudoabbreviata_G1113_S6032_1",
#              "Poa_lettermanii_G1230_S7434",
#              "Poa_laxa_G1268_Schonswetter8872_3",
#              "Poa_flexuosa_G0787_Brysting96117",
#              "Poa_compressa_G0209_6457",
#              "Poa_compressa_G2107_10338",
#              "AN18N065",
#              "Poa_hartzii_hartzii_G0246_6623_5",
#              "Poa_dolosa_G1244_S7495_1",
#              "Poa_diaphora_oxyglumis_G2103_10313",
#              "Poa_badensis_G1279_Hajkova2004_12",
#              "Poa_ligulata_G0903_JACA166095",
#              "Poa_flabellata_G1357_Wright9NSG",
#              "Poa_cookii_G1453_HennionGen1",
#              "Poa_supina_G1088_S5950_2",
#              "Poa_annua_G0164_6284" ,
#              "Poa_marcida_G2027_S5974",
#              "Poa_autumnalis_G1052_S4680_1" ,
#              "Poa_saltuensis_G0497_7043",
#              "Poa_wolfii_G1065_S5800"    
             )
#write.table(y, file = "~/poa/subsetted_tree_names.txt", sep = "\t", row.names = F, col.names = F)


# >> Convert tree to phylo  object ----
phy_tree <- as.phylo(final_tree) 

# >> Subset tree ----
keep_tree <- keep.tip(phy_tree, tip = species)

# >> Alter names ----
y <- as_tibble(keep_tree)

# Examples of edits:
# y$label[2] <- "Poa compressa Gillespie6457"
# y$label[3] <- "Poa compressa Gillespie10338"
# y$label[4] <- "Poa flexuosa Brysting96117"
# y$label[5] <- "Poa laxa Schonswetter8872-3"
# y$label[6] <- "Poa lettermanii Soreng7434"
# y$label[7] <- "Poa pseudoabbreviata Soreng6032-1"

y_label_edits <- read.csv("~/poa/subsetted_tree_names.txt", sep = "\t", header = F)

y$label <- y_label_edits$V4

# convert back to tree
my_subtree <- as.phylo(y)

# >> Plot tree ----  
ggtree(my_subtree) +
  geom_treescale(x=0.18, y=0, fontsize = 2.5) +
  geom_tiplab(size = 3) + ggplot2::xlim(0, 0.25) +
  geom_hilight(node = 55, fill = "darkgray", alpha = 0.3, extend = 0.08) +
  geom_hilight(node = 61, fill = "darkgray", alpha = 0.3, extend = 0.09)

ggsave("~/poa/Trees/ITS.tree.09132021_subsettedtree.jpeg", width = 15, height = 20, dpi = 300, units = "cm" )

# ## alt versions --
# 
# # plain tree - just names
# ggtree(my_subtree) +
#   geom_treescale(x=0.18, y=0, fontsize = 2.5) +
#   geom_tiplab(size = 3) + ggplot2::xlim(0, 0.25)

ggsave("~/poa/Trees/ITS.tree.09132021_subsettedtree_plain.jpeg", width = 15, height = 20, dpi = 300, units = "cm" )
