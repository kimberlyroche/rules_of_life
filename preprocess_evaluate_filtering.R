source("include/R/general.R")

level <- "genus"
glom_data <- load_glommed_data(level=level, replicates=TRUE)
subsetted_data <- subset_samples(glom_data, sname %in% sname_list)

filter_data_alt <- function(data, count_threshold=100, sample_threshold=10, verbose=TRUE) {
  snames <- unique(sample_data(data)$sname)
  counts <- otu_table(data)@.Data
  total_counts <- sum(counts)
  # a taxon must appear in at least {count_threshold} abundance in at least {sample_threshold} samples per individual
  keep_indices <- rep(TRUE, ntaxa(data))
  for(sname in snames) {
    sname <<- sname
    indiv_data <- subset_samples(data, sname==sname)
    # these are the indices to remove!
    keep_indices_indiv <- apply(counts, 2, function(x) sum(x >= count_threshold) >= sample_threshold)
    keep_indices <- keep_indices & keep_indices_indiv
  }
  collapse_indices <- !keep_indices
  # collapse mitochondria too
  tt <- tax_table(data)@.Data
  collapse_indices[which(tt[,colnames(tt) == "family"] == "Mitochondria")] <- TRUE
  merged_data <- merge_taxa(data, which(collapse_indices == TRUE), 1)
  retained_counts <- sum(counts[,!collapse_indices])
  if(verbose) {
    cat("Collapsing",sum(collapse_indices == TRUE),"taxa of",ntaxa(data),"\n")
    collapse_tidx <- which(collapse_indices == TRUE)[1]
    cat("\tOther category collapsed into index:",collapse_tidx,sep=" ","\n")
    cat("\tTaxonomy of collapsed (BEFORE):",tax_table(merged_data)@.Data[collapse_tidx,],"\n")
    tax_table(merged_data)@.Data[collapse_tidx,] <- rep("Collapsed",7)
    cat("\tTaxonomy of collapsed (AFTER):",tax_table(merged_data)@.Data[collapse_tidx,],"\n")
    cat("\tCollapsed counts:",(total_counts-retained_counts),"of",total_counts,"(",(total_counts-retained_counts)/total_counts,"total )\n")
    cat("\tPercent zero-count in data set:",(1 - get_tiny_counts(merged_data, 1)),"\n")
  }
  return(merged_data)
}

# let's assay some filtering options
cts <- c(1, 3, 5, 10, 100)
sts <- c(2, 5, 10, 20) # these are percent of samples
for(count_threshold in cts) {
  for(sample_threshold in sts) {
    cat("Cutoff:",count_threshold,",",sample_threshold,"\n")
    test <- filter_data_alt(subsetted_data, count_threshold, sample_threshold, verbose=TRUE)
  }
}
