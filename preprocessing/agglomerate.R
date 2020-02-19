# this file collapses microbes taxonomically to the specified level using phyloseq

relative_path <- ".."

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  cat("Arguments: (level)")
  quit()
}

level <- args[1]

source(file.path(relative_path,"include/R/general.R"))

if(level != "species") {
  # start from species
  glom_data <- readRDS(file.path(relative_path,data_dir,"glom_data_species_reps_tree.rds"))
  glom_data <- tax_glom(glom_data, taxrank=level, NArm=FALSE)
  saveRDS(glom_data, file=file.path(relative_path,data_dir,paste("glom_data_",level,"_reps_tree.rds")))
}
