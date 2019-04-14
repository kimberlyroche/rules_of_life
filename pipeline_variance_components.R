# this is a wrapper that executes a variance components fitting for the 16S data

source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for task missing!\n")
}

task <- as.numeric(args[1])
# 1 : 16S, no noise
# 2 : 16S, noise matrix included
# 3 : metagenomics, no noise
# 4 : metagenomics, noise included

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(filtered)

if(task == 1) {
  estimate_variance_components(filtered, metadata, optim_it=1, use_individuals=50, include_residual=F)
}
if(task == 2) {
  estimate_variance_components(filtered, metadata, optim_it=1, use_individuals=50, include_residual=T)
}

if(task == 3 || task == 4) {
  # read in metagenomics/PiPhillin
  data.piphillin <- read_metagenomics(metadata)
  metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)

  if(task == 3) {
    estimate_variance_components(data.piphillin, metadata.metagenomics, optim_it=1, use_individuals=50, include_residual=F)
  }
  if(task == 4) {
    estimate_variance_components(data.piphillin, metadata.metagenomics, optim_it=1, use_individuals=50, include_residual=T)
  }
}

