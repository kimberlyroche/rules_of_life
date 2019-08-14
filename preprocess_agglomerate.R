# collapse taxa via phyloseq

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  cat("Arguments: (level)")
  quit()
}

level <- args[1]

source("include/general.R")

if(level != "species") {
  load(paste0(data_dir,"glom_data_species_reps.RData"))
  glom_data <- glom_counts(glom_data, level=level)
  save(glom_data, file=paste0(data_dir,"glom_data_",level,"_reps.RData"))
}
