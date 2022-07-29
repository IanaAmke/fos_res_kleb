library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)
library(tidyverse)

# read blast datasets
Egli_blast <- read.delim(file = "fos_blast/Egli_MIC_gene_hits.tsv", header = TRUE, sep = "\t")
Huynh_blast <- read.delim(file = "fos_blast/Huynh2020_gene_hits.tsv", header = TRUE, sep = "\t")
Spark_blast <- read.delim(file = "fos_blast/spark_KpSc_gene_hits.tsv", header = TRUE, sep = "\t")
Sands_BARNARDS_blast <- read.delim(file = "fos_blast/sands2021_BARNARDS_gene_hits.tsv", header = TRUE, sep = "\t")
Sands_arthropods_blast <- read.delim(file = "fos_blast/sands2021_arthropods_gene_hits.tsv", header = TRUE, sep = "\t")

# delete sequencing name and add sample name to match antibiogram data
#Huynh_blast_1 <- read.delim(file = "genome_fos_blast/Huynh2020_gene_hits.tsv", header = TRUE, sep = "\t") %>%
  #separate(name, into = "name", sep = "[_]", extra = "drop")

#Huynh_blast_final <- inner_join(Huynh_blast_1, Huynh_sample_alias, by = c("name" = "run_accession"), 
                                #keep = FALSE) %>%
  #relocate(sample_alias, .after = name) %>%
  #select(sample_alias:length.to) %>%
  #rename(name = sample_alias)

#write_delim(Huynh_blast_final, file = "Huynh2020_gene_hits.tsv", delim = "\t")

# filter for 100% ID and query 
Egli_blast_100 <- Egli_blast %>%
  filter(Percent.ID == 100, Percent.query == 100)

Huynh_blast_100 <- Huynh_blast %>%
  filter(Percent.ID == 100, Percent.query == 100)

Sands_BARNARDS_blast_100 <- Sands_BARNARDS_blast %>%
  filter(Percent.ID == 100, Percent.query == 100)

Sands_arthropods_blast_100 <- Sands_arthropods_blast %>%
  filter(Percent.ID == 100, Percent.query == 100)

Spark_blast_100 <- Spark_blast %>%
  filter(Percent.ID == 100, Percent.query == 100)

# combine datasets with 100% ID and query vertically
fos_blast_100 <- rbind(Egli_blast_100, Huynh_blast_100, Spark_blast_100, Sands_BARNARDS_blast_100)

# load antibiogram data and combine the two datasets
antibiogram_data_df1 <- read.delim(file = "data/antibiogram_combined.tsv", header = TRUE, 
                                   sep = "\t", strip.white = TRUE)

antibiogram_data_df2 <- read.delim(file ="data/antibiogram_combined_2.tsv", header = TRUE,
                                   sep = "\t", strip.white = TRUE)

antibiogram_data <- full_join(antibiogram_data_df1, antibiogram_data_df2) 

# filter for fosfomycin in antibiogram data
# use latest AST guidelines (2022)
fos_antibiogram <- antibiogram_data %>%
  mutate(Antibiotic = case_when(Antibiotic == "Fosfomycin" ~ "fosfomycin",
                                TRUE ~ Antibiotic)) %>%
  filter(Antibiotic == "fosfomycin") %>%
  mutate(strain = Sample.Name) %>%
  relocate(strain, .after = Sample.Name) %>%
  separate(strain, into = "strain", sep = "[.]", extra = "drop") 

# merge fos_blast_100 and fos_antibiogram datasets
fosgene_fosfomycin <- inner_join(fos_antibiogram, fos_blast_100, by = c("Sample.Name" = "name"), 
                                keep = FALSE)

# bar plots showing total number of samples with fos genes
fosgene_fosfomycin %>%
  ggplot(aes(index, fill = gene)) + geom_bar(position = position_dodge(0.7), width = 0.5) +  
  labs(title = "Fosfomycin Resistance Genes Distribution",x = "Dataset", 
       y = "Count", fill = "Fosfomycin resistance genes")

# UPSET PLOTS
# tidy data
tidy_fosgene_fosfomycin <- fosgene_fosfomycin %>% 
  as_tibble(rownames = NA)

# MIC broth violin plots for individual res patterns
# filter for MIC broth
micbroth_tidy_fosgene_fosfomycin <- tidy_fosgene_fosfomycin %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

micbroth_tidy_fosgene_fosfomycin %>%
  group_by(gene) %>% 
  summarize(gene_list = list(gene), Measurement = Measurement) %>%
  ggplot(aes(gene_list, factor(Measurement))) + geom_violin() + geom_count(aes(color = gene)) + 
  labs(x = "fos genes", y = "Fosfomycin MIC broth measurement (mg/L)", 
       color = "fos genes", size = "Number of strains") + scale_x_upset()

# Disk diffusion violin plots for individual res patterns
# filter for disk diffusion
diskdiff_tidy_fosgene_fosfomycin <- tidy_fosgene_fosfomycin %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

diskdiff_tidy_fosgene_fosfomycin %>%
  group_by(gene) %>% 
  summarize(gene_list = list(gene), Measurement = Measurement) %>%
  ggplot(aes(gene_list, Measurement)) + geom_violin() + geom_count(aes(color = gene)) + 
  labs(x = "fos genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       color = "fos genes", size = "Number of strains") + scale_x_upset()

# MIC agar violin plots for individual res patterns
# filter for MIC agar
#micagar_tidy_fosgene_fosfomycin <- tidy_fosgene_fosfomycin %>%
  #filter(Laboratory.Typing.Method == "MIC (agar dilution)")

#micagar_tidy_fosgene_fosfomycin %>%
  #group_by(gene) %>% 
  #summarize(gene_list = list(gene), Measurement = Measurement) %>%
  #ggplot(aes(gene_list, Measurement)) + geom_violin() + geom_count(aes(color = gene)) + 
  #labs(x = "fos genes", y = "Fosfomycin MIC agar measurement (mg/L)", 
       #color = "fos genes", size = "Number of strains") + scale_x_upset()


