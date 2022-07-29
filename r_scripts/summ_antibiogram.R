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
  ggplot(aes(x = index, fill = Laboratory.Typing.Method)) +
  geom_bar(alpha = 0.8) +
  labs(title = "Total sample count per dataset", y = "Total sample count", fill = "Lab typing method") +
  scale_fill_brewer(palette = "Accent") + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# summary of testing standards used
summ_test_std <- fos_antibiogram %>%
  group_by(Testing.standard.year.or.version, index, Laboratory.Typing.Method) %>%
  count(Testing.standard) %>%
  mutate(MIC.breakpoint = case_when(Testing.standard.year.or.version == "v.10.0" ~ "S <= 322; R > 322",
                                    Testing.standard.year.or.version == "" ~ "",
                                    TRUE  ~ "S <= 32; R > 32"),
         Diskdiff.breakpoint = case_when(Testing.standard.year.or.version == "2015" ~ "In preparation",
                                         Testing.standard.year.or.version == "" ~ "",
                                         TRUE  ~ "S >= 24; R < 24")) %>%
  relocate(Testing.standard.year.or.version, MIC.breakpoint, Diskdiff.breakpoint, .after = Testing.standard) %>%
  rename(Dataset = index, Total = n)

# summary plots of total sample counts vs measurements for different lab typing methods  
# Disk diffusion
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  ggplot(aes(Measurement, fill = index)) + geom_bar() + 
  labs(title = "Disk Diffusion Measurements", x = "Measurement (mm)", y = "Frequency", fill = "Dataset") 

# summary of MIC measurements
# MIC agar dilution 
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  ggplot(aes(x = factor(Measurement), fill = index)) + geom_bar() + 
  labs(title = "MIC Agar Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Dataset") 

# MIC broth dilution 
fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  ggplot(aes(x = factor(Measurement), fill = index)) + geom_bar() + 
  labs(title = "MIC Broth Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Dataset") + 
  scale_x_discrete(limits = factor(2^seq(2, 9, by = 1)))


# plot lab testing methods for one dataset against the other datasets

# create colour variable 
cols <- c("other_dataset" = "purple", "current_dataset" = "cyan1")

# disk diffusion
dataset_disk_diff <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") 

dataset_disk_diff_noindex <- dataset_disk_diff %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_disk_diff, aes(Measurement)) + geom_bar(data = dataset_disk_diff_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin Disk Diffusion", x = "Measurement (mm)", fill = "Dataset")

# MIC (agar dilution)
dataset_MIC_agar <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)")

dataset_MIC_agar_noindex <- dataset_MIC_agar %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_agar, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_agar_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin MIC Agar Dilution", x = "Measurement (mg/L)", fill = "Dataset") 

# MIC (broth dilution)
dataset_MIC_broth <- fos_antibiogram %>%
  group_by(index) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)")

dataset_MIC_broth_noindex <- dataset_MIC_broth %>%
  ungroup() %>%
  select(Sample.Name:Laboratory.Typing.Method.Version.or.Reagent)

ggplot(dataset_MIC_broth, aes(factor(Measurement))) + geom_bar(data = dataset_MIC_broth_noindex, aes(fill = "other_dataset"), alpha = 0.5) +
  geom_bar(aes(fill = "current_dataset")) + facet_wrap(~ index, ncol = 1) + 
  labs(title = "Fosfomycin MIC Broth Dilution", x = "Measurement (mg/L)", fill = "Dataset") +
  scale_x_discrete(limits = factor(2^seq(3, 9, by = 1)))
