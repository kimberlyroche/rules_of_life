source("include/R/GP.R")

print_what <- list("genus"=c(1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 15, 17, 37))
print_what <- list("genus"=c(28, 36, 95))

for(level in names(print_what)) {
  idx <- print_what[[level]]
  fn <- paste0(data_dir,level,"_filtered.rds")
  if(!file.exists(fn)) {
    glom_data <- load_glommed_data(level=level, replicates=TRUE)
    subsetted_data <- subset_samples(glom_data, sname %in% over_50)
    data <- filter_data(subsetted_data, level=level, verbose=FALSE)
    saveRDS(data, fn)
  } else {
    data <- readRDS(fn)
  }
  cat("No. taxa:",ntaxa(data),"\n")
  tt <- tax_table(data)
  rownames(tt) <- NULL
  print(tt[idx,])
}

