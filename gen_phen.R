library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)
library(tidyverse)

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

write_delim(fos_antibiogram, file = "fos_antibiogram.tsv", delim = "\t")

# load kleborate data and combine the two datasets
kleborate_genotype_df1 <- read.delim(file = "data/kleb_results_concat.tsv", header = TRUE, 
                                 sep = "\t", strip.white = TRUE)

kleborate_genotype_df2 <- read.delim(file = "data/kleb_results_concat_2.tsv", header = TRUE,
                                     sep = "\t", strip.white = TRUE)
  
kleborate_genotype <- full_join(kleborate_genotype_df1, kleborate_genotype_df2)

# merge genotype and phenotype datasets
fos_ab_kleb <- full_join(fos_antibiogram, kleborate_genotype, keep = FALSE) %>%
  relocate(strain, species, species_match, .after = Sample.Name) %>%
  select(Sample.Name:O_locus_missing_genes, Fcyn_acquired, truncated_resistance_hits, 
         spurious_resistance_hits) %>%
  filter(Antibiotic == "fosfomycin") %>%
  mutate(Fcyn_acquired = case_when(Fcyn_acquired == "-" ~ "NKRD",
                                   TRUE ~ Fcyn_acquired),
         truncated_resistance_hits = case_when(truncated_resistance_hits == "-" ~ "NKRD",
                                               TRUE ~ truncated_resistance_hits),
         spurious_resistance_hits = case_when(spurious_resistance_hits == "-" ~ "NKRD",
                                              TRUE ~ spurious_resistance_hits)) 

write_delim(fos_ab_kleb, file = "fos_ab_kleb.tsv", delim = "\t")
 
# CREATE UPSET PLOTS
# convert tidy df (tibble)
tidy_fos_ab_kleb <- fos_ab_kleb %>% 
  as_tibble(rownames = NA) %>%
  select(Sample.Name:num_resistance_genes, Fcyn_acquired, truncated_resistance_hits, spurious_resistance_hits)

# create disk diffusion plots
# disk diffusion violin plot for all res patterns
# filter for disk diffusion
diskdiff_tidy_fos_ab_kleb_all <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(28, 29, 30)) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRD'))

diskdiff_tidy_fos_ab_kleb_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, Measurement)) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       size = "Number of strains") + scale_x_upset()  

# disk diffusion violin plots for individual res patterns
# filter for disk diffusion
diskdiff_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

# disk diffusion vs Fcyn_acquired 
diskdiff_tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       color = "Fosfomycin", size = "Number of strains") + scale_x_upset()

# disk diffusion vs truncated res hits
diskdiff_tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Trunc, Measurement)) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  labs(x = "Truncated resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       color = "Truncated res hits", size = "Number of strains") + scale_x_upset()
  
# disk diffusion vs spurious res hits
diskdiff_tidy_fos_ab_kleb %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, Measurement)) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  labs(x = "Spurious resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       color = "Spurious res hits", size = "Number of strains") + scale_x_upset()


# create MIC (agar dilution) plots
# MIC (agar dilution) violin plot for all res patterns
# filter for MIC (agar dilution)
micagar_tidy_fos_ab_kleb_all <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(28, 29, 30)) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRD'))

micagar_tidy_fos_ab_kleb_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, factor(Measurement))) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       size = "Number of strains") + scale_x_upset()

# MIC (agar dilution) violin plots for individual res patterns
# filter for MIC (agar dilution)
micagar_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  na.omit()

# MIC(agar dilution) vs Fcyn_acquired 
micagar_tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Fosfomycin", size = "Number of strains") +  scale_x_upset() 
  
# MIC (agar dilution) vs truncated res hits
micagar_tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Trunc, factor(Measurement))) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains") + scale_x_upset()

# MIC (agar dilution) vs spurious res hits
micagar_tidy_fos_ab_kleb %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains") + scale_x_upset()


# create MIC (broth dilution) plots
# MIC (broth dilution) violin plot for all res patterns
# filter for MIC (broth dilution)
micbroth_tidy_fos_ab_kleb_all <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(28, 29, 30)) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRD'))

micbroth_tidy_fos_ab_kleb_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, factor(Measurement))) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       size = "Number of strains") + scale_y_discrete(limits = factor(2^seq(2, 8, by = 1))) + 
  scale_x_upset()


# MIC (broth dilution) violin plots for individual res patterns
# filter for MIC (broth dilution)
micbroth_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  na.omit()

# MIC(broth dilution) vs Fcyn_acquired 
micbroth_tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Fosfomycin", size = "Number of strains") + scale_x_upset()

# MIC (broth dilution) vs truncated res hits
micbroth_tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Trunc, factor(Measurement))) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains") + 
  scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()

# MIC (broth dilution) vs spurious res hits
micbroth_tidy_fos_ab_kleb%>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains") + 
  scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()

# summary of res patterns and res pattern types
summ_res_patterns <- tidy_fos_ab_kleb %>%
  gather(Res.Pattern, Res.Pattern.Type, c(28, 29, 30)) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRD')) %>%
  group_by(index, Res.Pattern, Res.Pattern.Type, Laboratory.Typing.Method, Resistance.phenotype) %>%
  count()

# bar plots of fos res genes vs total count
# total number of strains for diff res patterns
# Fcyn_acquired
fcyn_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>%
  filter(str_detect(Fcyn_acquired, 'fos')) %>%
  mutate(Fcyn.present = case_when(Fcyn_acquired != "NKRD" ~ "fos"))

fcyn_tidy_fos_ab_kleb %>%
  ggplot(aes(Fcyn.present, fill = Fcyn_acquired)) + geom_bar(position = position_dodge(0.7), width = 0.5) +  
  labs(x = "Fosfomycin Resistance Genes", y = "Count", fill = "Fosfomycin resistance genes") +
  theme(axis.text.x = element_blank())

# truncated resistance hits
truncated_tidy_fos_ab_kleb  <- tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>%
  filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  mutate(trunc.hits.present = case_when(truncated_resistance_hits != "NKRD" ~ "truncated"))

truncated_tidy_fos_ab_kleb %>%
  ggplot(aes(trunc.hits.present, fill = truncated_resistance_hits)) + geom_bar(position = position_dodge(1.0)) +  
  labs(x = "Truncated Resistance Hits", y = "Count", fill = "Truncated resistance hits") +
  theme(axis.text.x = element_blank())

# spurious resistance hits
spurious_tidy_fos_ab_kleb  <- tidy_fos_ab_kleb %>%
  group_by(spurious_resistance_hits) %>%
  filter(str_detect(spurious_resistance_hits, 'fos')) %>%
  mutate(spurious.hits.present = case_when(spurious_resistance_hits != "NKRD" ~ "spurious"))

spurious_tidy_fos_ab_kleb %>%
  ggplot(aes(spurious.hits.present, fill = spurious_resistance_hits)) + geom_bar(position = position_dodge(1.0)) +  
  labs(x = "Spurious Resistance Hits", y = "Count", fill = "Spurious resistance hits") +
  theme(axis.text.x = element_blank())

# total number of strains for diff res patterns under diff lab typing methods
# create colours for resistant and susceptible strains
cols <- c("resistant" = "purple", "susceptible" = "cyan1")

# fcyn_acquired
fos_fcyn_present <- tidy_fos_ab_kleb %>%
  filter(str_detect(Fcyn_acquired, 'fos')) %>%
  mutate(Fcyn.present = case_when(Fcyn_acquired != "NKRD" ~ Fcyn_acquired)) %>%
  na.omit()

fos_fcyn_present %>%
  ggplot(aes(index, fill = Fcyn.present)) + geom_bar(position = "dodge", width = 0.25) + 
  facet_wrap(~ Laboratory.Typing.Method, nrow = 1) + labs(x = "Fosfomycin resistance genes", fill = "fos genes")

# truncated resistance hits
fos_trunc_present <- tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  mutate(Truncated.present = case_when(truncated_resistance_hits != "NKRD" ~ truncated_resistance_hits)) %>%
  na.omit()

fos_trunc_present %>%
  ggplot(aes(index, fill = Truncated.present)) + geom_bar(position = "dodge", width = 0.5) + 
  facet_wrap(~ Laboratory.Typing.Method, nrow = 1) + labs(x = "Truncated resistance hits", fill = "Truncated resistance hits")

# spurious resistance hits
fos_spur_present <- tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos')) %>%
  mutate(Spurious.present = case_when(spurious_resistance_hits != "NKRD" ~ spurious_resistance_hits)) %>%
  na.omit()

fos_spur_present %>%
  ggplot(aes(index, fill = Spurious.present)) + geom_bar(position = "dodge") + 
  facet_wrap(~ Laboratory.Typing.Method, nrow = 1) + labs(x = "Spurious resistance hits", fill = "Spurious resistance hits") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))