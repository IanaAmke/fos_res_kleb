library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)
library(tidyverse)

# load antibiogram data and combine the two datasets
antibiogram_data_df1 <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                                   sep = "\t", strip.white = TRUE)

antibiogram_data_df2 <- read.delim(file ="antibiogram_combined (1).tsv", header = TRUE,
                                   sep = "\t", strip.white = TRUE)

antibiogram_data <- full_join(antibiogram_data_df1, antibiogram_data_df2) 

# filter for fosfomycin in antibiogram data
# use latest AST guidelines (2022)
fos_antibiogram <- antibiogram_data %>%
  filter(Antibiotic == "fosfomycin") %>%
  mutate(Resistance.phenotype = case_when(Measurement < 6 & Laboratory.Typing.Method == "Disk diffusion" ~ "invalid", 
                                          Measurement < 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "resistant",
                                          Measurement >= 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "susceptible",
                                          Measurement > 8 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "resistant",
                                          Measurement <= 8 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "susceptible",
                                          Measurement > 8 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "resistant",
                                          Measurement <= 8 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "susceptible"))

# load kleborate data and combine the two datasets
kleborate_genotype_df1 <- read.delim(file = "kleb_results_concat.tsv", header = TRUE, 
                                 sep = "\t", strip.white = TRUE)

kleborate_genotype_df2 <- read.delim(file = "kleb_phenotype_all.tsv", header = TRUE,
                                     sep = "\t", strip.white = TRUE)
  

kleborate_genotype <- full_join(kleborate_genotype_df1, kleborate_genotype_df2)

# merge genotype and phenotype datasets
fos_ab_kleb <- full_join(fos_antibiogram, kleborate_genotype, 
                         by = c("Sample.Name" = "strain"), keep = TRUE) %>%
  relocate(strain, species, species_match, .after = Sample.Name) %>%
  mutate(Fcyn_acquired = coalesce(Fcyn_acquired, fosfomycin_gene)) %>%
  select(Sample.Name:O_locus_missing_genes, Fcyn_acquired, truncated_resistance_hits, 
         spurious_resistance_hits) %>%
  mutate(Fcyn_acquired = case_when(Fcyn_acquired == "-" ~ "NKRD",
                                   TRUE ~ Fcyn_acquired),
         truncated_resistance_hits = case_when(truncated_resistance_hits == "-" ~ "NKRD",
                                               TRUE ~ truncated_resistance_hits),
         spurious_resistance_hits = case_when(spurious_resistance_hits == "-" ~ "NKRD",
                                              TRUE ~ spurious_resistance_hits)) %>%
  filter(str_detect(Fcyn_acquired, 'fos|NKRD'))
  

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
       linetype = "Testing std cutoff", color = "Fosfomycin", size = "Number of strains") + 
  scale_x_upset()  

# disk diffusion violin plots for individual res patterns
# filter for disk diffusion
diskdiff_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

# disk diffusion vs Fcyn_acquired 
diskdiff_tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Testing std cutoff", color = "Fosfomycin", size = "Number of strains") +
  # scale_color_hue(labels = c("Known resistance gene", "No known resistance gene")) + 
  scale_color_manual(values = c("brown", "purple"), 
                     labels = c("Known resistance gene", "No known resistance gene")) +
  scale_linetype_manual(values = "dotted") + scale_x_upset()

# disk diffusion vs truncated res hits
diskdiff_tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Trunc, Measurement)) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Truncated resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Testing std cutoff", color = "Truncated res hits", size = "Number of strains") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()
  
# disk diffusion vs spurious res hits
diskdiff_tidy_fos_ab_kleb %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, Measurement)) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Spurious resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Zone diameter breakpoint", color = "Spurious res hits", size = "Number of strains") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()


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
  filter(Laboratory.Typing.Method == "MIC (agar dilution)")

# MIC(agar dilution) vs Fcyn_acquired 
micagar_tidy_fos_ab_kleb %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  #geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Fosfomycin", size = "Number of strains") +  scale_x_upset() 
  
# MIC (agar dilution) vs truncated res hits
micagar_tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Trunc, factor(Measurement))) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  geom_hline(aes(yintercept = 3, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()

# MIC (agar dilution) vs spurious res hits
micagar_tidy_fos_ab_kleb %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 3, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()


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
       size = "Number of strains") + scale_y_discrete(limits = factor(2^seq(3, 8, by = 1))) + scale_x_upset()

# MIC (broth dilution) violin plots for individual res patterns
# filter for MIC (broth dilution)
micbroth_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

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
  geom_hline(aes(yintercept = 2, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()

# MIC (broth dilution) vs spurious res hits
micbroth_tidy_fos_ab_kleb%>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 2, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()

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
trunc_tidy_fos_ab_kleb  <- tidy_fos_ab_kleb %>%
  group_by(truncated_resistance_hits) %>%
  filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  mutate(trunc.hits.present = case_when(truncated_resistance_hits != "NKRD" ~ "truncated"))

trunc_tidy_fos_ab_kleb %>%
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

# disk diffusion
# fosfomycin resistance genes
diskdiff_tidy_fos_ab_kleb %>%
  filter(str_detect(Fcyn_acquired, 'fos')) %>%
  mutate(Fcyn.present = case_when(Fcyn_acquired != "NKRD" ~ Fcyn_acquired)) %>%
  ggplot(aes(Fcyn.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.3) + 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "Disk Diffusion", x = "Fosfomycin resistance genes")

# truncated resistance hits - no truncated resistance hits
#diskdiff_tidy_fos_ab_kleb %>%
  #filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  #mutate(Truncated.present = case_when(truncated_resistance_hits != "NKRD" ~ truncated_resistance_hits)) %>%
  #ggplot(aes(Truncated.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.5)+ 
  #scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  #labs(title = "Disk Diffusion", x = "Truncated resistance hits")

# spurious resistance hits
diskdiff_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos')) %>%
  mutate(Spurious.present = case_when(spurious_resistance_hits != "NKRD" ~ spurious_resistance_hits)) %>%
  ggplot(aes(Spurious.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.5) + 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "Disk Diffusion", x = "Spurious resistance hits")

# MIC (agar dilution)
# fosfomycin resistance genes - no Fcyn_acquired entries
#micagar_tidy_fos_ab_kleb %>%
  #filter(str_detect(Fcyn_acquired, 'fos')) %>%
  #mutate(Fcyn.present = case_when(Fcyn_acquired != "NKRD" ~ Fcyn_acquired)) %>%
  #ggplot(aes(Fcyn.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.5) + 
  #scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  #labs(title = "MIC agar", x = "Fosfomycin resistance genes")

# truncated resistance hits
micagar_tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  mutate(Truncated.present = case_when(truncated_resistance_hits != "NKRD" ~ truncated_resistance_hits)) %>%
  ggplot(aes(Truncated.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.3)+ 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "MIC agar", x = "Truncated resistance hits")

# spurious resistance hits
micagar_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos')) %>%
  mutate(Spurious.present = case_when(spurious_resistance_hits != "NKRD" ~ spurious_resistance_hits)) %>%
  ggplot(aes(Spurious.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.6) + 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "MIC agar", x = "Spurious resistance hits")

# MIC (broth dilution)
# fosfomycin resistance genes - no Fcyn_acquired entries
#micbroth_tidy_fos_ab_kleb %>%
  #filter(str_detect(Fcyn_acquired, 'fos')) %>%
  #mutate(Fcyn.present = case_when(Fcyn_acquired != "NKRD" ~ Fcyn_acquired)) %>%
  #ggplot(aes(Fcyn.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.5) + 
  #scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  #labs(title = "MIC agar", x = "Fosfomycin resistance genes")

# truncated resistance hits
micbroth_tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos')) %>%
  mutate(Truncated.present = case_when(truncated_resistance_hits != "NKRD" ~ truncated_resistance_hits)) %>%
  ggplot(aes(Truncated.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.3) + 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "MIC broth", x = "Truncated resistance hits")

# spurious resistance hits
micbroth_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos')) %>%
  mutate(Spurious.present = case_when(spurious_resistance_hits != "NKRD" ~ spurious_resistance_hits)) %>%
  ggplot(aes(Spurious.present, fill = Resistance.phenotype)) + geom_bar(position = "dodge", width = 0.5) + 
  scale_color_manual(values = cols, aesthetics = "fill", name = "Resistance phenotype") + 
  labs(title = "MIC broth", x = "Spurious resistance hits")

# total number of strains for w/ or w/o fos genes that are res or sus
# disk diffusion
# Fcyn_acquired
diskdiff_tidy_fos_ab_kleb %>%
  filter(str_detect(Fcyn_acquired, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(Fcyn_acquired == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      Fcyn_acquired == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = Fcyn_acquired)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "Disk Diffusion", x = "Genotype vs Phenotype")

# truncated resistance hits
diskdiff_tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(truncated_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      truncated_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = truncated_resistance_hits)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "Disk Diffusion", x = "Genotype vs Phenotype")

# spurious resistance hits
diskdiff_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(spurious_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      spurious_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = spurious_resistance_hits)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "Disk Diffusion", x = "Genotype vs Phenotype")


# MIC (agar)
# Fcyn_acquired
micagar_tidy_fos_ab_kleb %>%
  filter(str_detect(Fcyn_acquired, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(Fcyn_acquired == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      Fcyn_acquired == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = Fcyn_acquired)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "MIC (agar dilution)", x = "Genotype vs Phenotype")

# truncated resistance hits
micagar_tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(truncated_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      truncated_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = truncated_resistance_hits)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "MIC (agar dilution)", x = "Genotype vs Phenotype")

# spurious resistance hits
micagar_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(spurious_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      spurious_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = spurious_resistance_hits)) + geom_bar(position = "dodge", width = 0.7) + 
  labs(title = "MIC (agar dilution)", x = "Genotype vs Phenotype")


# MIC (broth)
# Fcyn_acquired
micbroth_tidy_fos_ab_kleb %>%
  filter(str_detect(Fcyn_acquired, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(Fcyn_acquired == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      Fcyn_acquired == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      Fcyn_acquired != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = Fcyn_acquired)) + geom_bar(position = "dodge", width = 0.3) + 
  labs(title = "MIC (broth dilution)", x = "Genotype vs Phenotype")

# truncated resistance hits
micbroth_tidy_fos_ab_kleb %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(truncated_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      truncated_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      truncated_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = truncated_resistance_hits)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "MIC (broth dilution)", x = "Genotype vs Phenotype")

# spurious resistance hits
micbroth_tidy_fos_ab_kleb %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRD')) %>%
  mutate(gen.phen.pattern = case_when(spurious_resistance_hits == "NKRD" & Resistance.phenotype == "resistant" ~ "NKRD+r",
                                      spurious_resistance_hits == "NKRD" & Resistance.phenotype == "susceptible" ~ "NKRD+s",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "resistant" ~ "fos+r",
                                      spurious_resistance_hits != "NKRD" & Resistance.phenotype == "susceptible" ~ "fos+s")) %>%
  ggplot(aes(gen.phen.pattern, fill = spurious_resistance_hits)) + geom_bar(position = "dodge", width = 0.5) + 
  labs(title = "MIC (broth dilution)", x = "Genotype vs Phenotype")
