library(here)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)

# read newick file containing tree data
gene_tree <- read.newick("extracted_seq/filtered_fosgenes_NR_seq_2_midroot.nwk")

# prepare tree metadata
# read unique hits file
fos_fcyn_uniquehits <- read.delim("data/fosgenes_antibiogram_uniqueHits.tsv", header = TRUE, sep = "\t")

# read blast file (duplicates queried against non-redundant sequences)
blast <- read.delim("extracted_seq/fos_genes_dups.tsv", header = TRUE, sep = "\t")

# filter for 100% ID
blast100 <- blast %>%
  filter(Percent.ID == 100) %>%
  rename(gene = Subject.seq.ID)
#separate(Query.seq.ID, into = "Query.seq.ID", sep = "[_]", extra = "drop")

# create gene column
fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  mutate(gene = paste(strain, contig, fos.gene, sep = "_")) %>%
  rename(gene.copy.number = contig.number)

# combine unique hits file and blast output to fill the gene column with extracted fos genes
fosgene_fosfomycin <- left_join(fos_fcyn_uniquehits, blast100, by = c("gene" = "Query.seq.ID")) %>%
  select(Sample.Name:gene.y) %>%
  group_by(Sample.Name) %>%
  mutate(gene = case_when(gene.y %in% NA ~ gene, 
                          str_extract(gene, 'fos.*$') == str_extract(gene.y, 'fos.*$') ~ gene.y,
                          fos.gene == "fosA|MK043329.1" & gene.y == "SB5418_2_fosA6|ARO_3004111" ~ gene.y)) %>%
  drop_na(gene) %>%
  select(- gene.y) %>%
  relocate(gene, .after = index) %>%
  group_by(Sample.Name) %>%
  mutate(gene.copy.number = n()) %>%
  relocate(gene.copy.number, .after = gene) %>%
  mutate(strain = paste(strain, contig, fos.gene, sep = "_")) %>%
  group_by(Sample.Name)

filtered_fosgene_fosfomycin <- fosgene_fosfomycin %>%
  select(strain, Laboratory.Typing.Method, Measurement, Measurement.units, gene, fos.gene, contig, percent.ID, gene.copy.number)

# check for missing entries in extracted sequences that are in the unique hits file
#diff_unique_extracted <- fos_fcyn_uniquehits %>%
  #anti_join(filtered_fosgene_fosfomycin, by = c("gene" = "strain"))

#write_delim(filtered_fos_fcyn_uniquehits, file = "fosgenes_copynum_pheno.tsv", delim = "\t")

# read tree metadata
fos_gene_metadata <- filtered_fosgene_fosfomycin %>%
  ungroup()

# data for gene copy number
df1 <- fos_gene_metadata %>%
  group_by(gene, gene.copy.number) %>%
  count() %>%
  rename(total = n)

df1$gene.copy.number <- as.factor(df1$gene.copy.number)

# data for phenotypes
# mic broth dilution
df2 <- fos_gene_metadata %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  group_by(gene) %>%
  count(Measurement) %>%
  rename(measurement = Measurement, count = n)

df2$measurement <- as.factor(df2$measurement)

# disk diffusion
df3 <- fos_gene_metadata %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  select(gene, Measurement) %>%
  group_by(gene) %>%
  count(Measurement) %>%
  rename(measurement = Measurement, count = n)

df3$measurement <- as.factor(df3$measurement)

# adjust order
df2$measurement <- factor(df2$measurement, 
                                levels=c("16", "32", "64", "128", "256"))

df3$measurement <- factor(df3$measurement, 
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
p0 <- ggtree(gene_tree) +
  geom_tiplab(align = TRUE, linesize = 0.1, size = 2.0) +
  geom_hilight(data=kleborate_genes, mapping=aes(node=node),
               extendto=0.9, alpha=0.2, fill="darkmagenta",
               to.bottom = TRUE) +
  geom_tree(color="black", size=0.6) +
  theme_tree() +
  theme(
    legend.box = "vertical",
    legend.position = "right"
  ) 

p0

# add tree metadata
p1 <- p0 + 
  geom_fruit(data=df1, geom=geom_tile,
             mapping=aes(y=gene, x=gene.copy.number, fill=as.integer(total), group=total),
             color = "white", offset = 0.04, size = 0.01, width=0.01,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5,
                              title = "Gene copy \nnumber", title.height = 0.005,
                              hjust = 0.05), pwidth = 0.015
  ) + 
  scale_fill_gradient(low = "#f5faf6", high = "#04b31b", name = "Genomes \nper gene copy") +
  new_scale_fill() +
  geom_fruit(data=df2, geom=geom_tile,
             mapping=aes(y=gene, x=measurement, fill=count),
             color = "white", offset = 0.01, size = 0.01, width=0.01,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "MIC broth dilution", title.height = 0.005,
                              hjust = 0.05), pwidth = 0.06
  ) + 
  scale_fill_gradient2(low = "#faebeb", high = "#f50a0a", name = "MIC broth dilution\ncount", na.value = "grey") +
  new_scale_fill() +
  geom_fruit(data=df3, geom=geom_tile,
             mapping=aes(y=gene, x=measurement, fill=count),
             color = "white", offset = 0.01, size = 0.01,
             axis.params=list(axis="x", text.angle = 90, text.size = 1.5, 
                              title = "Disk diffusion", title.height = 0.005,
                              hjust = 0.05), pwidth = 0.15
  ) +
  scale_fill_continuous(low = "#ebebfa", high = "#0a0af5", name = "Disk diffusion\ncount") +
  geom_treescale(fontsize=2, linesize=0.3, x=0, y=100) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm")
  ) +
  coord_cartesian(clip = "off")

# view tree
p1

# save tree to pdf
ggsave("fos_gene_tree.pdf", width = 50, height = 50, units = "in", limitsize = FALSE)
