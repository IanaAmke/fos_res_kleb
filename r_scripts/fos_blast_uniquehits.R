library(here)
library(ggupset)
library(tidyverse)
library(gghighlight)

# read unique hits data files
Egli_uniquehits <- read.delim("fos_blast_unique_hits/Egli_MIC_uniqueHits.tsv", header = TRUE, sep = "\t")
Huynh_uniquehits <- read.delim("fos_blast_unique_hits/Huynh2020_uniqueHits.tsv", header = TRUE, sep = "\t") 
SPARK_uniquehits <- read.delim("fos_blast_unique_hits/spark_KpSc_uniqueHits.tsv", header = TRUE, sep = "\t")
sands_arthropods_uniquehits <- read.delim("fos_blast_unique_hits/sands2021_arthropods_uniqueHits.tsv", header = TRUE, sep = "\t")
sands_BARNARDS_uniquehits <- read.delim("fos_blast_unique_hits/sands2021_BARNARDS_uniqueHits.tsv", header = TRUE, sep = "\t")

# extract the run accession number 
#Huynh_uniquehits <- read.delim("fos_blast_unique_hits/Huynh2020_uniqueHits.tsv", header = TRUE, sep = "\t") %>%
  #separate(strain, into = "strain", sep = "[_]", extra = "drop")

# read Huynh2020_sample_alias.tsv file 
#Huynh_samplealias <- read.delim("fos_blast_unique_hits/Huynh2020_sample_alias.tsv", header = TRUE, sep = "\t")

# merge unique hits file to sample alias file to extract strain name
#Huynh_uniquefinal <- inner_join(Huynh_uniquehits, Huynh_samplealias, by = c("strain" = "run_accession"), keep = FALSE) %>%
  #relocate(sample_alias, .after = strain) %>%
  #select(sample_alias:coverage) %>%
  #rename(strain = sample_alias)

# save updated dataframe with correct strain names to tsv file
#write_delim(Huynh_uniquefinal, file = "Huynh2020_uniqueHits.tsv", delim = "\t")

# combine unique hits datasets vertically
fos_uniquehits <- rbind(Egli_uniquehits, Huynh_uniquehits, sands_arthropods_uniquehits, sands_BARNARDS_uniquehits, SPARK_uniquehits)

# read fos_antibiogram data file
fos_antibiogram <- read.delim(file = "data/fos_antibiogram.tsv", header = TRUE, sep = "\t")

# merge fos_blast and fos_antibiogram datasets
fos_fcyn_uniquehits <- inner_join(fos_antibiogram, fos_uniquehits, by = c("Sample.Name" = "strain"), 
                                 keep = FALSE)

# create dataframe with contig.number column that counts number of fos genes in a genome
fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  group_by(strain) %>%
  add_count() %>%
  rename(contig.number = n) 

# relocate strain after index and contig.number after contig
fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  reduce(
    .x = list(c('strain','index'), c('contig.number','contig')), 
    .f = ~ relocate(.x, .y[1], .after = .y[2]),
    .init = fos_fcyn_uniquehits
  )

# check if fos genes are in the same contigs
fos_fcyn_uniquehits_samecontig <- fos_fcyn_uniquehits %>%
  group_by(strain, contig) %>%
  add_count() %>%
  ungroup()

# delete one of fos genes in same contig (with lower %ID)
fos_fcyn_uniquehits_samecontig_filtered <- subset(fos_fcyn_uniquehits_samecontig, !duplicated(fos_fcyn_uniquehits_samecontig[c('strain', 'contig')]))

# delete contig.number and n columns
fos_fcyn_uniquehits_samecontig_filtered <- subset(fos_fcyn_uniquehits_samecontig_filtered, select = -c(n))

# rename fos_fcyn_uniquehits_samecontig_filtered to fos_fcyn_uniquehits and count gene copies
fos_fcyn_uniquehits <- fos_fcyn_uniquehits_samecontig_filtered %>%
  group_by(strain) %>%
  add_count() %>%
  mutate(contig.number = paste(n)) %>%
  select(- n)

# save fos_fcyn_uniquehits to tsv file
#write_delim(fos_fcyn_uniquehits, file = "fosgenes_antibiogram_uniqueHits.tsv", delim = "\t")

# bar plots showing total number of samples with fos genes
fos_fcyn_uniquehits %>%
  ggplot(aes(fos.gene)) + geom_bar(fill = "#0000CC") +  
  labs(title = "Fosfomycin Resistance Genes Distribution",x = "fos gene", 
       y = "Count", fill = "Fosfomycin resistance genes") + 
  # add labels on top of bars
  geom_text(aes(label = ..count..), stat = "count", nudge_y = 12) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# frequency table of fos gene distribution
fos_genes_freq <- fos_fcyn_uniquehits %>%
  group_by(fos.gene) %>%
  summarise(count = n()) 

# bar plots showing total number of samples with fos genes per dataset 
fos_fcyn_uniquehits %>%
  ggplot(aes(index, fill = fos.gene)) + geom_bar() +  
  labs(title = "Fosfomycin Resistance Genes Distribution Per Dataset",x = "Dataset", 
       y = "Count", fill = "Fosfomycin resistance genes") 

# bar plots showing total number of samples with fos genes per lab typing method 
# MIC broth dilution
fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  ggplot(aes(factor(Measurement), fill = fos.gene)) + geom_bar() +
  # ECOFF cut-off of 128 mg/L
  geom_vline(aes(xintercept = 4, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "Fosfomycin Resistance Genes Distribution for MIC (Broth Dilution)",x = "MIC (broth dilution) measurement (mg/L)", 
       y = "Count", fill = "Fosfomycin resistance genes", linetype = "MIC ECOFF cut-off") +
  scale_linetype_manual(values = "dotted") 

# MIC agar dilution
#fos_fcyn_uniquehits %>%
  #filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  #ggplot(aes(factor(Measurement), fill = fos.gene)) + geom_bar() +  
  # ECOFF cut-off of 128 mg/L
  #geom_vline(aes(xintercept = 7, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  #labs(title = "Fosfomycin Resistance Genes Distribution for MIC (Agar Dilution)",x = "MIC (agar dilution) measurement (mg/L)", 
       #y = "Count", fill = "Fosfomycin resistance genes", linetype = "MIC ECOFF cut-off") +
  #scale_linetype_manual(values = "dotted")

# Disk diffusion
fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  ggplot(aes(Measurement, fill = fos.gene)) + geom_bar() + 
  # ECOFF cut-off of 15 mm
  geom_vline(aes(xintercept = 15, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "Fosfomycin Resistance Genes Distribution for Disk Diffusion",x = "Disk diffusion measurement (mm)", 
       y = "Count", fill = "Fosfomycin resistance genes", linetype = "Disk diffusion \nECOFF cut-off") +
  scale_linetype_manual(values = "dotted")

# bar plots showing total number of samples with fos genes per lab typing method and dataset
# MIC broth dilution
fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  ggplot(aes(factor(Measurement), fill = fos.gene)) + geom_bar() +  
  labs(title = "Fosfomycin Resistance Genes Distribution Per Dataset for MIC (Broth Dilution)",x = "MIC (broth dilution) measurement (mg/L)", 
       y = "Count", fill = "Fosfomycin resistance genes") + 
  facet_wrap(~ index, nrow = 1)

# MIC agar dilution
#fos_fcyn_uniquehits %>%
  #filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  #ggplot(aes(factor(Measurement), fill = fos.gene)) + geom_bar() +  
  #labs(title = "Fosfomycin Resistance Genes Distribution Per Dataset for MIC (Agar Dilution)",x = "MIC (agar dilution) measurement (mg/L)", 
       #y = "Count", fill = "Fosfomycin resistance genes") + 
  #facet_wrap(~ index, nrow = 1)

# Disk diffusion
fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  ggplot(aes(factor(Measurement), fill = fos.gene)) + geom_bar() +  
  labs(title = "Fosfomycin Resistance Genes Distribution Per Dataset for Disk Diffusion",x = "Disk diffusion measurement (mm)", 
       y = "Count", fill = "Fosfomycin resistance genes") + 
  facet_wrap(~ index, nrow = 1)

# UPSET PLOTS
# MIC broth violin plots for individual res patterns
# ECOFF cut-offs added to infer phenotype
# filter for MIC broth
micbroth_fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

micbroth_fos_fcyn_uniquehits %>%
  group_by(Sample.Name, Measurement) %>% 
  summarize(gene_list = list(fos.gene)) %>%
  ggplot(aes(gene_list, Measurement)) + geom_violin() + geom_count() + 
  # ECOFF cut-off of 128 mg/L
  geom_hline(aes(yintercept = 128, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "Fos Gene Pattern for MIC (Broth Dilution)", x = "fos genes", y = "Fosfomycin MIC broth measurement (mg/L)",
       size = "Number of strains", linetype = "MIC ECOFF cut-off") +
  scale_linetype_manual(values = "dotted") +
  scale_x_upset() +
  theme_combmatrix(combmatrix.panel.line.size = 0.8)

# MIC agar violin plots for individual res patterns
# filter for MIC agar
#micagar_fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  #filter(Laboratory.Typing.Method == "MIC (agar dilution)")

#micagar_fos_fcyn_uniquehits %>%
  #group_by(Sample.Name, Measurement) %>% 
  #summarize(gene_list = list(fos.gene)) %>%
  #ggplot(aes(gene_list, Measurement)) + geom_violin() + geom_count() + 
  # ECOFF cut-off of 128 mg/L
  #geom_hline(aes(yintercept = 128, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  #labs(title = "Fos Gene Pattern for MIC (Agar Dilution)", x = "fos genes", y = "Fosfomycin MIC agar measurement (mg/L)",
  #size = "Number of strains", linetype = "MIC ECOFF cut-off") +
  #scale_linetype_manual(values = "dotted") +
  #scale_x_upset() +
  #theme_combmatrix(combmatrix.panel.line.size = 0.8)

# Disk diffusion violin plots for individual res patterns
# filter for disk diffusion
diskdiff_fos_fcyn_uniquehits <- fos_fcyn_uniquehits %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

diskdiff_fos_fcyn_uniquehits %>%
  group_by(Sample.Name, Measurement) %>% 
  summarize(gene_list = list(fos.gene)) %>%
  ggplot(aes(gene_list, Measurement)) + geom_violin() + geom_count() + 
  # ECOFF cut-off of 15 mm
  geom_hline(aes(yintercept = 15, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "Fos Gene Pattern for Disk Diffusion", x = "fos genes", y = "Fosfomycin disk diffusion measurement (mm)",
       size = "Number of strains", linetype = "Disk diffusion \nECOFF cut-off") +
  scale_linetype_manual(values = "dotted") +
  scale_x_upset() +
  theme_combmatrix(combmatrix.panel.line.size = 0.8)
  