library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)

# read newick file containing tree data
gene_tree <- read.newick("phylogenetic_trees/fos_genes_NR.nwk")

# read tree metadata
fos_gene_metadata <- read.delim("data/fos_genes_metadata.tsv", header = TRUE, sep = "\t")

# data for phenotypes
# mic broth dilution
df1 <- fos_gene_metadata %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  select(gene, Measurement) %>%
  group_by(gene) %>%
  count(Measurement) %>%
  rename(measurement = Measurement, count = n)

df1$measurement <- as.factor(df1$measurement)

# disk diffusion
df2 <- fos_gene_metadata %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  select(gene, Measurement) %>%
  group_by(gene) %>%
  count(Measurement) %>%
  rename(measurement = Measurement, count = n)

df2$measurement <- as.factor(df2$measurement)

# adjust order
df1$measurement <- factor(df1$measurement, 
                                levels=c("16", "32", "64", "128", "256"))

df2$measurement <- factor(df2$measurement, 
                                 levels=c("6", "13", "14", "15", "16",
                                          "17", "18", "19", "20", "21",
                                          "22", "23", "24", "25", "27"))

# highlight nodes with fos genes from kleborate
kleborate_genes <- as.tibble(gene_tree) %>%
  left_join(
    fos_gene_metadata %>% select(gene), by = c("label" = "gene")
  ) %>%
  mutate(gene = paste(label)) %>%
  filter(grepl('[*]$', gene)) %>%
  unique()
 
# plot tree and highlight strains with genes from kleborate
p <- ggtree(gene_tree) +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 2.0) +
  geom_hilight(data=kleborate_genes, mapping=aes(node=node),
               extendto=0.97, alpha=0.2, fill="darkmagenta",
               to.bottom = TRUE) +
  geom_tree(color="black", size=0.6) +
  theme_tree() +
  theme(
    legend.box = "vertical",
    legend.position = "right"
  ) 

p 

# add tree metadata
p <- p + 
  geom_fruit(data=df1, geom=geom_tile,
             mapping=aes(y=gene, x=measurement, fill=count),
             color = "white", offset = 0.04, size = 0.02,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "MIC broth dilution", title.height = 0.005,
                              hjust = 0.05), pwidth = 0.1
  ) + 
  scale_fill_continuous(low = "#faebeb", high = "#f50a0a", name = "MIC broth dilution\ncount") +
  new_scale_fill() +
  geom_fruit(data=df2, geom=geom_tile,
             mapping=aes(y=gene, x=measurement, fill=count),
             color = "white", offset = 0.06, size = 0.02 ,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "Disk diffusion", title.height = 0.005,
                              hjust = 0.05), pwidth = 0.15
  ) +
  scale_fill_continuous(low = "#ebebfa", high = "#0a0af5", name = "Disk diffusion\ncount") +
  geom_treescale(fontsize=2, linesize=0.3, x=0, y=100) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  ) +
  coord_cartesian(clip = "off")

# view tree
p

# save tree to pdf
ggsave("fos_gene_tree.pdf", width = 50, height = 50, units = "in", limitsize = FALSE)
