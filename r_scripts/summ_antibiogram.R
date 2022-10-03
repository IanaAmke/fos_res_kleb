library(here)
library(dplyr)
library(ggplot2)

# load antibiogram data and combine the two datasets
antibiogram_data_df1 <- read.delim(file = "data/antibiogram_combined.tsv", header = TRUE, 
                               sep = "\t", strip.white = TRUE)

antibiogram_data_df2 <- read.delim(file ="data/antibiogram_combined_2.tsv", header = TRUE,
                                   sep = "\t", strip.white = TRUE)

antibiogram_data <- full_join(antibiogram_data_df1, antibiogram_data_df2) 

# filter for fosfomycin 
# use latest AST guidelines (2022)
fos_antibiogram <- antibiogram_data %>%
  mutate(Antibiotic = case_when(Antibiotic == "Fosfomycin" ~ "fosfomycin",
                                TRUE ~ Antibiotic)) %>%
  filter(Antibiotic == "fosfomycin") 

# summary of samples within each dataset
summ_dataset <- fos_antibiogram %>%
  group_by(index, Laboratory.Typing.Method, Laboratory.Typing.Platform) %>%
  count() %>%
  rename(Dataset = index, Total = n)

# summary plot of total sample count by dataset
fos_antibiogram %>%
  group_by(index) %>%
  ggplot(aes(x = index)) +
  geom_bar(fill = "#ff8000", alpha = 0.6) +
  geom_text(aes(label = ..count..), stat = "count", nudge_y = 60) +
  labs(title = "Total sample count per dataset", y = "Total sample count", fill = "Lab typing method") +
  scale_fill_brewer(palette = "Accent") + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# summary plots of total sample counts vs measurements for different lab typing methods  
# ECOFF cut-offs used due to there being no K. pneumoniae breakpoints for fosfomycin
# Disk diffusion
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  ggplot(aes(Measurement, fill = index)) + geom_bar() + 
  # ECOFF cut-off of 15 mm
  geom_vline(aes(xintercept = 15, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "Disk Diffusion Measurements", x = "Measurement (mm)", y = "Frequency", fill = "Dataset", linetype = "Disk diffusion \nECOFF cut-off") +
  scale_linetype_manual(values = "dotted")

# summary of MIC measurements
# MIC agar dilution 
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  ggplot(aes(x = factor(Measurement), fill = index)) + geom_bar() + 
  # ECOFF cut-off of 128 mg/L
  geom_vline(aes(xintercept = 7, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "MIC Agar Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Dataset", linetype = "MIC ECOFF cut-off") +
  scale_linetype_manual(values = "dotted")

# MIC broth dilution 
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  ggplot(aes(x = factor(Measurement), fill = index)) + geom_bar() + 
  # ECOFF cut-off of 128 mg/L
  geom_vline(aes(xintercept = 6, linetype = "EUCAST"), alpha = 0.6, color = "black") +
  labs(title = "MIC Broth Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Dataset", linetype = "MIC ECOFF cut-off") +
  scale_linetype_manual(values = "dotted") +
  scale_x_discrete(limits = factor(2^seq(2, 9, by = 1)))


# plot lab testing methods for one dataset against the other datasets

# create colour variable 
cols <- c("other dataset" = "purple", "current dataset" = "cyan1")

# disk diffusion
dataset_disk_diff <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") 

dataset_disk_diff_noindex <- dataset_disk_diff %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_disk_diff, aes(Measurement)) + geom_bar(data = dataset_disk_diff_noindex, aes(fill = "other dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin Disk Diffusion Measurements per Dataset", x = "Measurement (mm)", fill = "Dataset")

# MIC (agar dilution)
dataset_MIC_agar <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)")

dataset_MIC_agar_noindex <- dataset_MIC_agar %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_agar, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_agar_noindex, aes(fill = "other dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin MIC Agar Dilution Measurements per Dataset", x = "Measurement (mg/L)", fill = "Dataset") 

# MIC (broth dilution)
dataset_MIC_broth <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

dataset_MIC_broth_noindex <- dataset_MIC_broth %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_broth, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_broth_noindex, aes(fill = "other dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin MIC Broth Dilution Measurements per Dataset", x = "Measurement (mg/L)", fill = "Dataset") +
  scale_x_discrete(limits = factor(2^seq(3, 9, by = 1)))
