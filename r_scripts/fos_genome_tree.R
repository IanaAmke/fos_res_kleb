library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)
#library(wesanderson)
library(colorspace)

# read newick file containing tree data
tree <- read.newick("phylogenetic_trees/fos_gene_with_dups.nwk")

# prepare tree metadata
fos_fcyn_uniquehits <- read.delim("data/fos_fcyn_uniquehits.tsv", header = TRUE, sep = "\t")

blast <- read.delim("data/fos_genes_dups.tsv", header = TRUE, sep = "\t")

blast100 <- blast %>%
  filter(Percent.ID == 100) %>%
  rename(gene = Subject.seq.ID)%>%
  separate(Query.seq.ID, into = "Query.seq.ID", sep = "[_]", extra = "drop")

fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  mutate(gene = paste(strain, contig, fos.gene, sep = "_"))

fosgene_fosfomycin <- full_join(fos_fcyn_uniquehits, blast100, by = c("Sample.Name" = "Query.seq.ID"), 
                                keep = FALSE) %>%
  select(Sample.Name:gene.y) %>%
  mutate(gene.y = ifelse(gene.y %in% NA, gene.x, gene.y)) %>%
  select(- gene.x) %>%
  rename(gene = gene.y) %>%
  relocate(gene, .after = index) %>%
  group_by(Sample.Name) %>%
  mutate(gene.copy.number = n()) %>%
  relocate(gene.copy.number, .after = gene) %>%
  mutate(strain = paste(strain, contig, fos.gene, sep = "_")) %>%
  group_by(Sample.Name)

filtered_fosgene_fosfomycin <- fosgene_fosfomycin %>%
  select(strain, Laboratory.Typing.Method, Measurement, Measurement.units, gene, gene.copy.number)

#write_delim(filtered_fosgene_fosfomycin, file = "fosgenes_copynum_pheno.tsv", delim = "\t")

# read tree metadata
fos_metadata <- read.delim("data/fosgenes_copynum_pheno.tsv", header = TRUE, sep = "\t")

# data for genes 
dat1 <- fos_metadata %>%
  select(strain, gene)

# data for gene copy number
dat2 <- fos_metadata %>%
  select(strain, gene.copy.number) %>%
  mutate(total.gene.number = paste(gene.copy.number))

dat2$total.gene.number <- as.factor(dat2$total.gene.number)

# data for phenotypes
# mic broth dilution
dat3 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  select(strain, Measurement) %>%
  mutate(MIC.measurements = paste(Measurement))

dat3$MIC.measurements <- as.factor(dat3$MIC.measurements)

# disk diffusion
dat4 <- fos_metadata %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  select(strain, Measurement) %>%
  mutate(D.df.measurements = paste(Measurement))

dat4$D.df.measurements <- as.factor(dat4$D.df.measurements)

# adjust order
dat2$total.gene.number <- factor(dat2$total.gene.number, 
                     levels=c("1", "2", "3", "4"))

dat3$MIC.measurements <- factor(dat3$MIC.measurements, 
                     levels=c("16", "32", "64", "128", "256"))

dat4$D.df.measurements <- factor(dat4$D.df.measurements, 
                                levels=c("6", "13", "14", "15", "16",
                                         "17", "18", "19", "20", "21",
                                         "22", "23", "24", "25", "27"))

# highlight nodes with instrisic fosA5 gene
kleborate_genes <- as.tibble(tree) %>%
  left_join(
    fos_metadata %>% select(strain, gene), by = c("label" = "strain")
  ) %>%
  filter(grepl('[*]$', gene))

# plot tree and highlight strains with intrinsic fosA5
p0 <- ggtree(tree) %<+% dat2 +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 0.5) +
  geom_hilight(data=kleborate_genes, mapping=aes(node=node),
               extendto=0.933, alpha=0.3, fill="darkorange",
               size=0.05, to.bottom = TRUE, align="right")

p0

# add tree metadata
# colour palettes
#colour_palette <- wes_palette("Moonrise3", type = "continuous")
#colour_palette_2 <- wes_palette("BottleRocket2", type = "continuous")
#colour_palette_3 <- wes_palette("Zissou1", 100, type = "continuous")

p1 <- p0 + 
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=gene.copy.number),
             color = "white", offset = 0.05, size = 0.02, width = 0.02,
             axis.params=list(axis="x", text.size = 0, title = "Gene copy number", title.height = 0.01,
                              line.alpha = 0), pwidth = 0.1
             ) + 
  scale_fill_continuous(low = "#e4e5f7", high = "#0a12f5", name = "Gene copy\nnumber") +
  new_scale_fill() +
  geom_fruit(data=dat3, geom=geom_tile,
             mapping=aes(y=strain, x=MIC.measurements, fill=Measurement),
             color = "white", offset = 0.03, size = 0.02 ,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "MIC broth dilution", title.height = 0.01,
                              hjust = 0.05), pwidth = 0.15
             ) +
  scale_fill_continuous(low = "#e4f7e9", high = "#03fc45", name = "MIC broth dilution\nmeasurements") +
  new_scale_fill() +
  geom_fruit(data=dat4, geom=geom_tile,
             mapping=aes(y=strain, x=D.df.measurements, fill=Measurement),
             color = "white", offset = 0.02, size = 0.02 ,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "Disk diffusion", title.height = 0.01,
                              hjust = 0.05), pwidth = 0.4
  ) +
  scale_fill_continuous(low = "#f7e4e4", high = "#fa0505", name = "Disk diffusion\nmeasurements") +
  geom_treescale(fontsize=2, linesize=0.3, x=0, y=1000) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm")
  ) +
  coord_cartesian(clip = "off")

# view tree
p1

# save tree to pdf
ggsave("fos_genome_tree.pdf", width = 50, height = 50, units = "in", limitsize = FALSE)
