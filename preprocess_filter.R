# collapse taxa via phyloseq

source("include/R/general.R")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  cat("Arguments: (level)")
  quit()
}

level <- args[1]

cat(paste0("Filtering data at the ",level," level...\n"))

count_threshold <- 5
sample_threshold <- 20

glom_data <- readRDS(paste0(data_dir,"glom_data_",level,"_reps_tree.rds"))
subsetted_data <- subset_samples(glom_data, sname %in% sname_list)
filtered_data <- filter_data(subsetted_data, level=level,
                             count_threshold=count_threshold, sample_threshold=sample_threshold,
                             verbose=TRUE)
