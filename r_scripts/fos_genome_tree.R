library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)

# read newick file containing tree data
tree <- read.newick("phylogenetic_trees/filtered_fosgenes_NR_ID_midpoint.nwk")

# prepare tree metadata
fos_fcyn_uniquehits <- read.delim("data/fosgenes_antibiogram_uniqueHits.tsv", header = TRUE, sep = "\t")








fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  mutate(gene = paste(strain, contig, fos.gene, sep = "_")) %>%
  rename(gene.copy.number = contig.number)



filtered_fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  select(strain, Laboratory.Typing.Method, Measurement, Measurement.units, gene, fos.gene, contig, percent.ID, gene.copy.number)

#write_delim(filtered_fos_fcyn_uniquehits, file = "fosgenes_copynum_pheno.tsv", delim = "\t")

# read tree metadata
fos_metadata <- filtered_fos_fcyn_uniquehits

# data for gene copy number
dat1 <- fos_metadata %>%
  select(gene, gene.copy.number) %>%
  mutate(total.gene.number = paste(gene.copy.number))

dat1$total.gene.number <- as.factor(dat1$total.gene.number)

# data for phenotypes
# mic broth dilution
dat2 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  select(gene, Measurement) %>%
  mutate(MIC.measurements = paste(Measurement))

dat2$MIC.measurements <- as.factor(dat2$MIC.measurements)

# disk diffusion
dat3 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  select(gene, Measurement) %>%
  mutate(D.df.measurements = paste(Measurement))

dat3$D.df.measurements <- as.factor(dat3$D.df.measurements)

# adjust order
dat1$total.gene.number <- factor(dat1$total.gene.number, 
                     levels=c("1", "2", "3", "4"))

dat2$MIC.measurements <- factor(dat2$MIC.measurements, 
                     levels=c("16", "32", "64", "128", "256"))

dat3$D.df.measurements <- factor(dat3$D.df.measurements, 
                                levels=c("6", "13", "14", "15", "16",
                                         "17", "18", "19", "20", "21",
                                         "22", "23", "24", "25", "27"))

# highlight nodes with instrisic fosA5 gene
kleborate_genes <- as.tibble(tree) %>%
  left_join(
    fos_metadata %>% select(gene, fos.gene), by = c("label" = "gene")
  ) %>%
  filter(grepl('[*]$', fos.gene))

# plot tree and highlight strains with intrinsic fosA5
p0 <- ggtree(tree) %<+% dat1 +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 0.5) +
  geom_hilight(data=kleborate_genes, mapping=aes(node=node),
               extendto=0.933, alpha=0.3, fill="darkorange",
               size=0.05, to.bottom = TRUE, align="right")

p0

# add tree metadata


p1 <- p0 + 
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=gene.copy.number),
             color = "white", offset = 0.05, size = 0.02, width = 0.02,
             axis.params=list(axis="x", text.size = 0, title = "Gene copy number", title.height = 0.01,
                              line.alpha = 0), pwidth = 0.1
             ) + 
  scale_fill_continuous(low = "#ffffff", high = "#0a12f5", name = "Gene copy\nnumber", breaks = c(0,1,2),
                        limits = c(0,2)) +
  new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=gene, x=MIC.measurements, fill=Measurement),
             color = "white", offset = 0.03, size = 0.02 ,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "MIC broth dilution", title.height = 0.01,
                              hjust = 0.05), pwidth = 0.15
             ) +
  scale_fill_continuous(low = "#e4f7e9", high = "#03fc45", name = "MIC broth dilution\nmeasurements") +
  new_scale_fill() +
  geom_fruit(data=dat3, geom=geom_tile,
             mapping=aes(y=gene, x=D.df.measurements, fill=Measurement),
             color = "white", offset = 0.02, size = 0.02 ,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "Disk diffusion", title.height = 0.01,
                              hjust = 0.05), pwidth = 0.4
  ) +
  scale_fill_continuous(low = "#faebeb", high = "#f50a0a", name = "Disk diffusion\nmeasurements") +
  geom_treescale(fontsize=2, linesize=0.3, x=0, y=1000) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.03, "cm")
  ) +
  coord_cartesian(clip = "off")

# view tree
p1

# save tree to pdf
ggsave("fos_genome_tree.pdf", width = 50, height = 50, units = "in", limitsize = FALSE)
