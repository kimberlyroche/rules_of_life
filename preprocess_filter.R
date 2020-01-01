# collapse taxa via phyloseq

source("include/R/general.R")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  cat("Arguments: (level) (count threshold) (sample percent threshold)")
  quit()
}

level <- args[1]
count_threshold <- as.numeric(args[2])
sample_threshold <- as.numeric(args[3])

cat(paste0("Filtering data at the ",level," level...\n"))

glom_data <- readRDS(paste0(data_dir,"glom_data_",level,"_reps_tree.rds"))
subsetted_data <- subset_samples(glom_data, sname %in% sname_list)
filtered_data <- filter_data(subsetted_data, level=level,
                             count_threshold=count_threshold, sample_threshold=sample_threshold,
                             verbose=TRUE)
