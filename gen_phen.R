library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)

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

# create upset plot
# convert df into a tibble
tidy_fos_ab_kleb <- fos_ab_kleb %>% 
  as_tibble(rownames = NA)

# create disk diffusion plot

#disk_diff_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  #group_by(Laboratory.Typing.Method) %>%
  #filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  #gather(Resistance.Pattern, Resistance.Pattern.Type, c(52, 54, 74, 75)) %>%
  #relocate(Resistance.Pattern, Resistance.Pattern.Type, .after = num_resistance_genes) %>%
  #select(Sample.Name:Resistance.Pattern.Type) 

# filter for disk diffusion
disk_diff_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")

# convert "-" to no known res genes (NKRG)
disk_diff_tidy_fos_ab_kleb[disk_diff_tidy_fos_ab_kleb == "-"] <- "NKRG"

#res_patterns <- disk_diff_tidy_fos_ab_kleb %>%
  #select.list(multiple = TRUE)

#res_patterns_list <- res_patterns %>%
  #as.list()

# disk diffusion plot for fosfomycin resistance
ggplot(disk_diff_tidy_fos_ab_kleb, aes(Fcyn_acquired, Measurement)) +
  geom_violin() + geom_count(aes(color = Fcyn_acquired)) + 
  geom_hline(aes(yintercept = 24, linetype = "EUCAST"), alpha = 0.6, color = "black") + 
  labs(x = "Fosfomycin acquired genes", y = "Fosfomycin disk diffusion measurement (mm)", 
       linetype = "Testing std cutoff", color = "Fosfomycin", size = "Number of strains") +
  scale_color_hue(labels = c("Known resistance gene", "No known resistance gene")) + scale_color_manual(values = c("brown", "purple")) +
  scale_linetype_manual(values = "dotted")



