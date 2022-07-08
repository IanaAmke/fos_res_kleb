library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)

# load antibiogram data
antibiogram_data <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                               sep = "\t")

# filter for fosfomycin in antibiogram data
fos_antibiogram <- antibiogram_data %>%
  filter(Antibiotic == "fosfomycin") %>%
  mutate(Resistance.phenotype = case_when(Measurement < 6 & Laboratory.Typing.Method == "Disk diffusion" ~ "invalid", 
                                          Measurement < 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "resistant",
                                          Measurement >= 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "susceptible",
                                          Measurement > 84 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "resistant",
                                          Measurement <= 84 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "susceptible",
                                          Measurement > 84 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "resistant",
                                          Measurement <= 84 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "susceptible"))

# load kleborate data
kleborate_genotype <- read.delim(file = "kleb_results_concat.tsv", header = TRUE, 
                                 sep = "\t")

# merge genotype and phenotype datasets
fos_ab_kleb <- inner_join(fos_antibiogram, kleborate_genotype, 
                          by = c("Sample.Name" = "strain"), keep = TRUE) %>%
  relocate(strain, species, species_match, .after = Sample.Name)


# CREATE UPSET PLOTS

# convert df_1 into tidy_df_1 (tibble)
# create new df of res patterns and types
tidy_fos_ab_kleb_df1 <- fos_ab_kleb %>% 
  as_tibble(rownames = NA) %>%
  gather(Res.Pattern, Res.Pattern.Type, c(54, 74, 75)) %>%
  relocate(Res.Pattern, Res.Pattern.Type, .after = num_resistance_genes) %>%
  select(Sample.Name:Res.Pattern.Type) %>%
  group_by(Res.Pattern.Type) %>%
  summarize(Res.Patterns = list(Res.Pattern)) 

# convert "-" to no known res genes (NKRG)
tidy_fos_ab_kleb_df1[tidy_fos_ab_kleb_df1 == "-"] <- "NKRG"

# total count of strains with resistance patterns
tidy_fos_ab_kleb_df1 %>%
  ggplot(aes(x = Res.Patterns)) + geom_bar() + scale_x_upset(order_by = "degree") +
  labs(x = "Resistance Patterns", y = "Count")

# convert df_2 into tidy_df_2 (tibble)
tidy_fos_ab_kleb_df2 <- fos_ab_kleb %>% 
  as_tibble(rownames = NA)

# convert "-" to no known res genes (NKRG)
tidy_fos_ab_kleb_df2[tidy_fos_ab_kleb_df2 == "-"] <- "NKRG"

# create disk diffusion plots
# disk diffusion violin plot for all res patterns
# filter for disk diffusion
diskdiff_tidy_fos_ab_kleb_df2_all <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(54, 74, 75)) %>%
  relocate(Res.Pattern, Res.Pattern.Type, .after = num_resistance_genes) %>%
  select(Sample.Name:Res.Pattern.Type) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRG'))

diskdiff_tidy_fos_ab_kleb_df2_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, Measurement)) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Testing std cutoff", color = "Fosfomycin", size = "Number of strains") + 
  scale_x_upset()  

# disk diffusion violin plots for individual res patterns
# filter for disk diffusion
diskdiff_tidy_fos_ab_kleb_df2 <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

# disk diffusion vs Fcyn_acquired 
diskdiff_tidy_fos_ab_kleb_df2 %>%
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
diskdiff_tidy_fos_ab_kleb_df2 %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Trunc, Measurement)) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Truncated resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Testing std cutoff", color = "Truncated res hits", size = "Number of strains") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()
  
# disk diffusion vs spurious res hits
diskdiff_tidy_fos_ab_kleb_df2 %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Spurious, Measurement)) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Spurious resistance hits", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Zone diameter breakpoint", color = "Spurious res hits", size = "Number of strains") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()


# create MIC (agar dilution) plots
# MIC (agar dilution) violin plot for all res patterns
# filter for MIC (agar dilution)
micagar_tidy_fos_ab_kleb_df2_all <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(54, 74, 75)) %>%
  relocate(Res.Pattern, Res.Pattern.Type, .after = num_resistance_genes) %>%
  select(Sample.Name:Res.Pattern.Type) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRG'))

micagar_tidy_fos_ab_kleb_df2_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, factor(Measurement))) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       size = "Number of strains") + scale_x_upset()

# MIC (agar dilution) violin plots for individual res patterns
# filter for MIC (agar dilution)
micagar_tidy_fos_ab_kleb_df2 <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)")

# MIC(agar dilution) vs Fcyn_acquired 
micagar_tidy_fos_ab_kleb_df2 %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  #geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Fosfomycin", size = "Number of strains") +  scale_x_upset() 
  
# MIC (agar dilution) vs truncated res hits
micagar_tidy_fos_ab_kleb_df2 %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Trunc, factor(Measurement))) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  geom_hline(aes(yintercept = 3, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()

# MIC (agar dilution) vs spurious res hits
micagar_tidy_fos_ab_kleb_df2 %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 3, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (agar) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_x_upset()


# create MIC (broth dilution) plots
# MIC (broth dilution) violin plot for all res patterns
# filter for MIC (broth dilution)
micbroth_tidy_fos_ab_kleb_df2_all <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  gather(Res.Pattern, Res.Pattern.Type, c(54, 74, 75)) %>%
  relocate(Res.Pattern, Res.Pattern.Type, .after = num_resistance_genes) %>%
  select(Sample.Name:Res.Pattern.Type) %>%
  filter(str_detect(Res.Pattern.Type, 'fos|NKRG'))

micbroth_tidy_fos_ab_kleb_df2_all %>%
  group_by(Res.Pattern.Type) %>% 
  summarize(Res.Patterns = list(Res.Pattern), Measurement = Measurement) %>%
  ggplot(aes(Res.Patterns, factor(Measurement))) + geom_violin() + geom_count() + 
  labs(x = "Resistance genes", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       size = "Number of strains") + scale_y_discrete(limits = factor(2^seq(3, 8, by = 1))) + scale_x_upset()

# MIC (broth dilution) violin plots for individual res patterns
# filter for MIC (broth dilution)
micbroth_tidy_fos_ab_kleb_df2 <- tidy_fos_ab_kleb_df2 %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

# MIC(broth dilution) vs Fcyn_acquired 
micbroth_tidy_fos_ab_kleb_df2 %>%
  group_by(Fcyn_acquired) %>% 
  summarize(Fcyn = list(Fcyn_acquired), Measurement = Measurement) %>%
  ggplot(aes(Fcyn, Measurement)) + geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Fosfomycin", size = "Number of strains") + scale_x_upset()

# MIC (broth dilution) vs truncated res hits
micbroth_tidy_fos_ab_kleb_df2 %>%
  group_by(truncated_resistance_hits) %>% 
  summarize(Trunc = list(truncated_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(truncated_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Trunc, factor(Measurement))) + geom_violin() + geom_count(aes(color = truncated_resistance_hits)) + 
  geom_hline(aes(yintercept = 2, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Truncated resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Truncated res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()

# MIC (broth dilution) vs spurious res hits
micbroth_tidy_fos_ab_kleb_df2 %>%
  group_by(spurious_resistance_hits) %>% 
  summarize(Spurious = list(spurious_resistance_hits), Measurement = Measurement) %>%
  filter(str_detect(spurious_resistance_hits, 'fos|NKRG')) %>%
  ggplot(aes(Spurious, factor(Measurement))) + geom_violin() + geom_count(aes(color = spurious_resistance_hits)) + 
  geom_hline(aes(yintercept = 2, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(x = "Spurious resistance hits", y = "Fosfomycin MIC (broth) measurement (mg/L)", 
       color = "Spurious res hits", size = "Number of strains", linetype = "MIC breakpoint") + 
  scale_linetype_manual(values = "dotted") + scale_y_discrete(limits = factor(2^seq(2, 9, by = 1))) + scale_x_upset()
