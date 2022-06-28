library(here)
library(dplyr)
library(ggplot2)

# load antibiogram data
antibiogram_data <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                               sep = "\t")

# filter for fosfomycin 
fos_antibiogram <- antibiogram_data %>%
  filter(Antibiotic == "fosfomycin")

# summary of samples within each dataset
summ_dataset <- fos_antibiogram %>%
  group_by(index) %>%
  count() %>%
  rename(Dataset = index, Total = n)

# summary of testing standards used
summ_test_std <- fos_antibiogram %>%
  count(Testing.standard)
  
# filter for disk diffusion measurements
summ_Disk_diff <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  mutate(Resistance.phenotype = case_when(Measurement < 6 ~ "invalid", 
                                          Measurement < 24 ~ "resistant",
                                          Measurement >= 24 ~ "susceptible"))
# plot for disk diffusion frequency
ggplot(summ_Disk_diff, aes(Measurement, fill = Resistance.phenotype)) + geom_histogram(binwidth = 1, alpha = 0.6) + 
  labs(title = "Disk Diffusion Measurements", x = "Measurement (mm)", y = "Frequency", fill = "Phenotype") 

# summary of MIC measurements
# filter for MIC agar dilution measurements
summ_MIC_agar <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (agar dilution)") %>%
  mutate(Resistance.phenotype = case_when(Measurement > 84 ~ "resistant",
                                          Measurement <= 84 ~ "susceptible")) 
# plot for MIC agar dilution frequency
ggplot(summ_MIC_agar, aes(Measurement, fill = Resistance.phenotype)) + geom_histogram(binwidth = 80, alpha = 0.8) + 
  scale_x_continuous(breaks = seq(0, 2200, 200)) + scale_y_continuous(breaks = seq(0, 300, 50)) +
  labs(title = "MIC Agar Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Phenotype") 

# filter for MIC broth dilution measurements
summ_MIC_broth <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "MIC (broth dilution)") %>%
  mutate(Resistance.phenotype = case_when(Measurement > 84 ~ "resistant",
                                          Measurement <= 84 ~ "susceptible"))
# plot for MIC broth dilution frequency
ggplot(summ_MIC_broth, aes(Measurement, fill = Resistance.phenotype)) + geom_histogram(binwidth = 100, alpha = 0.8) + 
  labs(title = "MIC Broth Dilution Measurements", x = "Measurement (mg/L)", y = "Frequency", fill = "Phenotype")

# fosfomycin data grouped by index
summ_measurement_dataset <- fos_antibiogram %>%
  group_by(index) %>%
  mutate(Resistance.phenotype = case_when(Measurement < 6 & Laboratory.Typing.Method == "Disk diffusion" ~ "invalid", 
                                          Measurement < 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "resistant",
                                          Measurement >= 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "susceptible",
                                          Measurement > 84 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "resistant",
                                          Measurement <= 84 & Laboratory.Typing.Method == "MIC (agar dilution)" ~ "susceptible",
                                          Measurement > 84 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "resistant",
                                          Measurement <= 84 & Laboratory.Typing.Method == "MIC (broth dilution)" ~ "susceptible"))
# plot of phenotype by dataset
ggplot(summ_measurement_dataset, aes(x = index, fill = Laboratory.Typing.Method)) +
  geom_bar(alpha = 0.8) + scale_y_continuous(breaks = seq(0, 300, 50)) + facet_wrap(~ Resistance.phenotype) +
  labs(title = "Phenotype Frequency Per Dataset", x = "Dataset", y = "Total count", fill = "Lab typing method") +
  scale_fill_brewer(palette = "Accent")

  
