# this file filters low-abundance taxa into an "Other" category and saves the results back into
# the data directory

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  cat("Arguments: (level) (count threshold) (sample percent threshold)")
  quit()
}

level <- args[1]
count_threshold <- as.numeric(args[2])
sample_threshold <- as.numeric(args[3])

cat(paste0("Filtering data at the ",level," level...\n"))

if(level == "ASV" | is.null(level)) {
	glom_data <- read_data(replicates=FALSE)
} else {
	glom_data <- readRDS(file.path(relative_path,data_dir,paste0("glom_data_",level,"_reps_tree.rds")))
}
subsetted_data <- subset_samples(glom_data, sname %in% sname_list)
filtered_data <- filter_data(subsetted_data, level=level,
                             count_threshold=count_threshold, sample_threshold=sample_threshold,
                             verbose=TRUE)
