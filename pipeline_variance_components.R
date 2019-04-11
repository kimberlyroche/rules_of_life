# this is a wrapper that executes a variance components fitting for the 16S data

source("include.R")

glom_data <- load("glom_data_species.RData")
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
estimate_variance_components(filtered=filtered, optim_it=1)
