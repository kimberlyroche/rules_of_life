# this is a wrapper that executes a variance components fitting for the 16S data

source("include.R")

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(filtered)

estimate_variance_components(filtered, metadata, optim_it=1)

# read in metagenomics/PiPhillin
data.piphillin <- read_metagenomics(metadata)
metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata.metagenomics)

estimate_variance_components(data.piphillin, metadata.metagenomics, optim_it=1)
