library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)
library(wesanderson)

# read newick file containing tree data
tree <- read.newick("phylogenetic_trees/fos_genes_genome_lvl.nwk")

# read tree metadata
fos_metadata <- read.delim("data/fosgenes_copynum_pheno.tsv", header = TRUE, sep = "\t")

# data for genes 
dat1 <- fos_metadata %>%
  select(strain, gene)

# data for gene copy number
dat2 <- fos_metadata %>%
  select(strain, gene.copy.number) 

# data for phenotypes
# mic broth dilution
dat3 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  select(strain, Measurement) 

# count fos genes under each measurement 
#x <- dat3 %>%
#group_by(strain, Measurement) %>%
#count() %>%
#rename(count = n) %>%
#pivot_wider(names_from = Measurement, values_from = count, names_sort = TRUE) 

# replace all NA values with 0
#x[is.na(x)] <- 0

#dat3$MIC.measurements <- as.numeric(dat3$MIC.measurements)

# disk diffusion
dat4 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  select(strain, Measurement) 

#dat4$D.df.measurements <- as.numeric(dat4$D.df.measurements)

# mic broth dataset
# read file containing lab typing method measurements for heatmap
mic_measurements <- read.delim("data/micbroth_measurements.tsv", sep = "\t", row.names = 1)
colnames(mic_measurements) <- sub("^X", "", colnames(mic_measurements))
row_names_mic <- rownames(mic_measurements)
mic_measurements <- as.data.frame(sapply(mic_measurements, as.numeric))
rownames(mic_measurements) <- row_names_mic

# disk diffusion dataset
# read file containing lab typing method measurements for heatmap
diskdiff_measurements <- read.delim("data/diskdiff_measurements.tsv", sep = "\t", row.names = 1)
colnames(diskdiff_measurements) <- sub("^X", "", colnames(diskdiff_measurements))
row_names_diskdiff <- rownames(diskdiff_measurements)
diskdiff_measurements <- as.data.frame(sapply(diskdiff_measurements, as.numeric))
rownames(diskdiff_measurements) <- row_names_diskdiff

# gene copy number
gene_copy <- read.delim("data/fos_gene_copy_no.tsv", sep = "\t", row.names = 1)
colnames(gene_copy) <- sub("^X", "", colnames(gene_copy))
row_names_gene_copy <- rownames(gene_copy)
gene_copy <- as.data.frame(sapply(gene_copy, as.numeric))
rownames(gene_copy) <- row_names_gene_copy

# build tree 
p0 <- ggtree(tree) +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 0.5)

# build heatmap
# build heatmap
colour_palette <- wes_palette("Moonrise3", 100, type = "continuous")
colour_palette_2 <- wes_palette("GrandBudapest2", 100, type = "continuous")
colour_palette_3 <- wes_palette("Zissou1", 100, type = "continuous")

p1 <- gheatmap(p0, gene_copy,offset = 0.1, colnames_position = "bottom",
               colnames_offset_y = 0, hjust = 0, font.size = 2) +
  scale_fill_gradientn(colours = colour_palette, name = "Gene copy\nnumber")

p2 <- p1 + new_scale_fill()

p3 <- gheatmap(p2, diskdiff_measurements, offset = 1.15, colnames_position = "bottom",
               colnames_offset_y = 0, hjust = 0, font.size = 2) +
  scale_fill_gradientn(colours = colour_palette_2, name = "Disk diffusion\ncount")

p4 <- p3 + new_scale_fill()

p5 <- gheatmap(p4, mic_measurements, offset = 2.05, colnames_position = "bottom",
               colnames_offset_y = 0, hjust = 0, 
               font.size = 2) +
  scale_fill_gradientn(colours = colour_palette_3, name = "MIC broth dilution\ncount")

# view tree
p5

# save tree to pdf file
ggsave("fos_tree.pdf", width = 50, height = 50, units = "in", limitsize = FALSE)
