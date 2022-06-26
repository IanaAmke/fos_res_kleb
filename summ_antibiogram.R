library(here)
library(dplyr)
library(ggplot2)

# load antibiogram data
antibiogram_data <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                               sep = "\t")

# filter for fosfomycin 
fos_antibiogram <- antibiogram_data %>%
  filter(Antibiotic == "fosfomycin")

# total count for each of the lab typing methods
lab_typing <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  count() %>%
  rename(Lab_typing = Laboratory.Typing.Method, Total = n)

ggplot(lab_typing, aes(Lab_typing, Total, fill = Lab_typing)) +
  geom_col(alpha = 0.7) + labs(x = "Lab typing method", fill = "Lab Typing")
  
# summary for disk diffusion
summ_disk_diff <- fos_antibiogram %>%
  group_by(Laboratory.Typing.Method) %>%
  filter(Laboratory.Typing.Method == "Disk diffusion") %>%
  transmute(Sample.Name, Antibiotic, Measurement, Measurement.units, Laboratory.Typing.Method, 
            Invalid = Measurement < 6, Resistant = Measurement < 24, Susceptible = Measurement >= 24, 
            Testing.standard)


