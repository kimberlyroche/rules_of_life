# collapse taxa via phyloseq

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  cat("Arguments: (level)")
  quit()
}

level <- args[1]

source("include/R/general.R")

if(level != "species") {
  # old file format:
  #   load(paste0(data_dir,"glom_data_species_reps_tree.RData"))
  # this function seems to have been deleted in some distant version of this repo
  #   glom_data <- glom_counts(glom_data, level=level) # this function seems to have been deleted in some distant version of this repo
  glom_data <- readRDS(paste0(data_dir,"glom_data_species_reps_tree.rds"))
  glom_data <- tax_glom(glom_data, taxrank=level, NArm=FALSE)
  saveRDS(glom_data, file=paste0(data_dir,"glom_data_",level,"_reps_tree.rds"))
}
