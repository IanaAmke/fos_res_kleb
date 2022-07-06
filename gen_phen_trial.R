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

# create upset plot
# convert df into a tibble
tidy_fos_ab_kleb <- fos_ab_kleb %>% 
  as_tibble(rownames = NA) %>%
  gather(Res.Pattern, Res.Pattern.Type, c(52, 54, 74, 75)) %>%
  relocate(Res.Pattern, Res.Pattern.Type, .after = num_resistance_genes) %>%
  select(Sample.Name:Res.Pattern.Type)

tidy_fos_ab_kleb[tidy_fos_ab_kleb == "-"] <- NA

tidy_fos_ab_kleb_list <- tidy_fos_ab_kleb %>%
  group_by(Res.Pattern.Type) %>%
  summarize(Res.Patterns = list(Res.Pattern)) 

# total count of strains with resistance patterns
tidy_fos_ab_kleb_list %>%
  ggplot(aes(x = Res.Patterns)) + geom_bar() + scale_x_upset(order_by = "degree") +
  labs(x = "Resistance Patterns", y = "Count")

# disk diffusion + res genes
diskdiff_tidy_fos_ab_kleb <- tidy_fos_ab_kleb %>%
  filter(Laboratory.Typing.Method == "Disk diffusion")
  
  


