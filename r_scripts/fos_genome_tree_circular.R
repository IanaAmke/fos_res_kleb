library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)
#library(wesanderson)

# read newick file containing tree data
tree <- read.newick("phylogenetic_trees/filtered_fosgenes_NR_ID_midpoint.nwk")

# prepare tree metadata
fos_fcyn_uniquehits <- read.delim("data/fosgenes_antibiogram_uniqueHits.tsv", header = TRUE, sep = "\t")

fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  mutate(gene = paste(strain, contig, fos.gene, sep = "_")) %>%
  rename(gene.copy.number = contig.number)

filtered_fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  select(strain, Laboratory.Typing.Method, Measurement, Measurement.units, gene, fos.gene, contig, percent.ID, gene.copy.number)

# read tree metadata
fos_metadata <- filtered_fos_fcyn_uniquehits

# data for gene copy number
dat1 <- fos_metadata %>%
  select(gene, gene.copy.number) %>%
  mutate(kleborate = case_when(grepl('[*]$', gene) == TRUE ~ gene))

# highlight nodes with instrisic fosA5 gene
kleborate_genes <- as.tibble(tree) %>%
  left_join(
    fos_metadata %>% select(gene, fos.gene), by = c("label" = "gene")
  ) %>%
  filter(grepl('[*]$', fos.gene))

# plot tree and highlight strains with intrinsic fosA5
p0 <- ggtree(tree, layout = "circular") %<+% dat1 +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 0.5) +
  geom_hilight(data=kleborate_genes, mapping=aes(node=node),
               extendto=0.94, alpha=0.3, fill="darkorange",
               size=0.05, to.bottom = TRUE, align="right")

p0

p1 <- p0 + 
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=gene.copy.number),
             color = "white", offset = 0.08, size = 0.02, width = 0.02,
             axis.params=list(axis="x", text.size = 0, title.height = 0.01,
                              line.alpha = 0), pwidth = 0.1
  ) + 
  scale_fill_continuous(low = "#ffffff", high = "#0a12f5", name = "Gene copy\nnumber", breaks = c(0,1,2),
                        limits = c(0,2)) +
  geom_treescale(fontsize=2, linesize=0.3) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=16),
        legend.text=element_text(size=10),
        legend.spacing.y = unit(0.06, "cm")
  )

# view tree
p1

# save tree to pdf
ggsave("fos_genome_tree_circular.pdf", width = 20, height = 20, units = "in", limitsize = FALSE)

