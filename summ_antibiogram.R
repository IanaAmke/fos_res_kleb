library(here)
library(dplyr)
library(ggplot2)

# load antibiogram data
antibiogram_data <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, sep = "\t")

# filter for fosfomycin 
fos_antibiogram <- antibiogram_data %>%
  filter(Antibiotic == "fosfomycin")


  
