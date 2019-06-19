source("include.R")

load("glom_data_genus_reps.RData")
glom_data_family <- glom_counts(glom_data, level="family")
save(glom_data_family, file="glom_data_family_reps.RData")
