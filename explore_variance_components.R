# this is a wrapper that executes a variance components fitting for the 16S data

source("include/R/general.R")

include_residual <- FALSE

level <- "species"
data <- load_and_filter(level)
metadata <- read_metadata(data)

estimate_variance_components(filtered, metadata, optim_it=1, use_individuals=25, include_residual=include_residual)

# note: this can be run with metagenomics/PiPhillin as:
# data.piphillin <- read_metagenomics(metadata)
# metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)
# estimate_variance_components(data.piphillin, metadata.metagenomics, optim_it=1, use_individuals=25, include_residual=F)
