mutate_res_phenotype <- function() {
  require(dplyr)
  
  dat <- read.delim(file = "antibiogram_combined.tsv", header = TRUE, 
                    sep = "\t")
  
  dat %>%
    filter(Antibiotic == "fosfomycin") %>%
    mutate(Resistance.phenotype = case_when(Measurement < 6 & Laboratory.Typing.Method == "Disk diffusion" ~ "invalid", 
                                            Measurement < 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "resistant",
                                            Measurement >= 24 & Laboratory.Typing.Method == "Disk diffusion" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2015" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2015" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2017" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2017" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2018" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "2018" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "v9.0" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "v9.0" ~ "susceptible",
                                            Measurement > 322 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "v.10.0" ~ "resistant",
                                            Measurement <= 322 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "v.10.0" ~ "susceptible",
                                            Measurement > 82 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "" ~ "resistant",
                                            Measurement <= 82 & Laboratory.Typing.Method == "MIC (agar dilution)" & Testing.standard.year.or.version == "" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2015" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2015" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2017" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2017" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2018" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2018" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "v9.0" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "v9.0" ~ "susceptible",
                                            Measurement > 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2019" ~ "resistant",
                                            Measurement <= 32 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "2019" ~ "susceptible",
                                            Measurement > 322 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "v.10.0" ~ "resistant",
                                            Measurement <= 322 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "v.10.0" ~ "susceptible",
                                            Measurement > 82 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "" ~ "resistant",
                                            Measurement <= 82 & Laboratory.Typing.Method == "MIC (broth dilution)" & Testing.standard.year.or.version == "" ~ "susceptible"))
}

