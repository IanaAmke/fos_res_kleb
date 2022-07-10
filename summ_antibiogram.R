library(here)
library(dplyr)
library(ggplot2)

# load antibiogram data and combine the two datasets
antibiogram_data_df1 <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                               sep = "\t", strip.white = TRUE)

antibiogram_data_df2 <- read.delim(file ="antibiogram_combined (1).tsv", header = TRUE,
                                   sep = "\t", strip.white = TRUE)

antibiogram_data <- inner_join(antibiogram_data_df1, antibiogram_data_df2) 

# filter for fosfomycin 
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

# summary of samples within each dataset
summ_dataset <- fos_antibiogram %>%
  group_by(index) %>%
  count() %>%
  rename(Dataset = index, Total = n)

# summary of testing standards used
summ_test_std <- fos_antibiogram %>%
  group_by(Testing.standard.year.or.version, index) %>%
  count(Testing.standard) %>%
  mutate(MIC.breakpoint = case_when(Testing.standard.year.or.version == "v.10.0" ~ "S <= 322; R > 322",
                                    Testing.standard.year.or.version == "" ~ "",
                                    TRUE  ~ "S <= 32; R > 32"),
         Diskdiff.breakpoint = case_when(Testing.standard.year.or.version == "2015" ~ "In preparation",
                                         Testing.standard.year.or.version == "" ~ "",
                                         TRUE  ~ "S >= 24; R < 24")) %>%
  relocate(Testing.standard.year.or.version, MIC.breakpoint, Diskdiff.breakpoint, .after = Testing.standard) %>%
  rename(Dataset = index, Total = n)
  
# filter for disk diffusion measurements
summ_Disk_diff <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") 
  
# plot for disk diffusion frequency
ggplot(summ_Disk_diff, aes(Measurement, fill = Resistance.phenotype)) + geom_bar() + 
  labs(title = "Disk Diffusion Measurements", x = "Measurement (mm)", y = "Frequency", fill = "Phenotype") 

# summary of MIC measurements
# filter for MIC agar dilution measurements
summ_MIC_agar <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") 
   
# plot for MIC agar dilution frequency
ggplot(summ_MIC_agar, aes(x = factor(Measurement), fill = Resistance.phenotype)) + geom_bar() + 
  labs(title = "MIC Agar Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Phenotype") 

# filter for MIC broth dilution measurements
summ_MIC_broth <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") 

# plot for MIC broth dilution frequency
ggplot(summ_MIC_broth, aes(x = factor(Measurement), fill = Resistance.phenotype)) + geom_bar() + 
  labs(title = "MIC Broth Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Phenotype") + 
  scale_x_discrete(limits = factor(2^seq(2, 9, by = 1)))

# fosfomycin data grouped by index
summ_measurement_dataset <- fos_antibiogram %>%
  group_by(index) 
  
# plot of phenotype by dataset
dataset_plot <- ggplot(summ_measurement_dataset, aes(x = index, fill = Laboratory.Typing.Method)) +
  geom_bar(alpha = 0.8) + scale_y_continuous(breaks = seq(0, 300, 50)) + facet_wrap(~ Resistance.phenotype) +
  labs(title = "Phenotype Frequency Per Dataset", y = "Total sample count", fill = "Lab typing method") +
  scale_fill_brewer(palette = "Accent") 

dataset_plot + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# plot lab testing methods for one dataset against the other datasets

# create colour variable 
cols <- c("other_dataset" = "purple", "current_dataset" = "cyan1")

# disk diffusion
dataset_disk_diff <- summ_measurement_dataset %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") 

dataset_disk_diff_noindex <- dataset_disk_diff %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_disk_diff, aes(Measurement)) + geom_bar(data = dataset_disk_diff_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + geom_vline(aes(xintercept = 24, linetype = "EUCAST"), 
                                                                                       alpha = 0.6, color = "black") +
  labs(title = "Fosfomycin Disk Diffusion", x = "Measurement (mm)", fill = "Dataset", linetype = "Zone diameter breakpoint") +
  scale_linetype_manual(values = "dotted")

# MIC (agar dilution)
dataset_MIC_agar <- summ_measurement_dataset %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)")

dataset_MIC_agar_noindex <- dataset_MIC_agar %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_agar, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_agar_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + geom_vline(aes(xintercept = 3, linetype = "EUCAST"), 
                                                                                       alpha = 0.6, color = "black") +
  labs(title = "Fosfomycin MIC Agar Dilution", x = "Measurement (mg/L)", fill = "Dataset", linetype = "MIC breakpoint") +
  scale_linetype_manual(values = "dotted")

# MIC (broth dilution)
dataset_MIC_broth <- summ_measurement_dataset %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

dataset_MIC_broth_noindex <- dataset_MIC_broth %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_broth, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_broth_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + geom_vline(aes(xintercept = 2, linetype = "EUCAST"), 
                                                                                       alpha = 0.6, color = "black") +
  labs(title = "Fosfomycin MIC Broth Dilution", x = "Measurement (mg/L)", fill = "Dataset", linetype = "MIC breakpoint") +
  scale_linetype_manual(values = "dotted") + scale_x_discrete(limits = factor(2^seq(2, 9, by = 1)))
