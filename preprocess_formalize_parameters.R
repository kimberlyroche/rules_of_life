# get empirical estimates for variance to use in GP kernels

source("include/R/general.R")
source("include/R/data_transform.R")
source("include/R/visualization.R")

levels <- c("phylum", "family", "genus")

for(level in levels) {
  # estimate noise in replicates
  # we want variance between samples!
  #data <- load_and_filter(level)
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, level=level, verbose=FALSE)
  subsetted_data <- subset_samples(data, sname %in% sname_list)
  cat("Data loaded...\n")
  md <- read_metadata(subsetted_data)
  sample_status <- md$sample_status
  cat("Percent replicates:",round(sum(sample_status == 2)/length(sample_status), 3),"\n")
  
  # apply log ratio transform
  cat("Number of taxa at this level:",ntaxa(subsetted_data),"\n")
  alr_ref <- pick_alr_ref(subsetted_data)
  cat("ALR ref index (median abundance):",alr_ref,"\n")
  log_ratios <- apply_alr(subsetted_data, d=alr_ref)
  unique_sids <- unique(md[md$sample_status == 2,]$sid)
  replicate_var <- c()
  for(usid_idx in 1:length(unique_sids)) {
    # omit sketchy replicates with super low total counts (< 5K)
    retain_sids <- md$sample_id[md$sid == unique_sids[usid_idx]]
    sample_counts <- otu_table(subsetted_data)@.Data[md$sid == unique_sids[usid_idx],]
    # omit where total counts are below 5K bounds
    total_counts <- c(apply(sample_counts, 1, sum))
    retain_sids <- retain_sids[total_counts >= 5000]
    if(length(retain_sids) >= 3) {
      # these are samples (rows) x taxa (columns)
      lr <- log_ratios[(rownames(log_ratios) %in% retain_sids),]
      lr <- scale(lr, scale=FALSE)
      # mean-centered the mean of each vector's norm should be an estimate of variance between replicates
      # (just from formula for variance)
      sample_var <- cov(t(lr))
      replicate_var <- c(replicate_var, diag(sample_var))
    }
  }
  
  indiv_var <- c()
  for(indiv in unique(md$sname)) {
    retain_sids <- md[md$sname == indiv,]$sample_id
    if(length(retain_sids) > 3) {
      lr <- log_ratios[(rownames(log_ratios) %in% retain_sids),]
      lr <- scale(lr, scale=FALSE)
      sample_var <- cov(t(lr))
      indiv_var <- c(indiv_var, diag(sample_var))
    }
  }
  
  mean_total_var <- mean(indiv_var)
  mean_rep_var <- mean(replicate_var)
  #mean_after_reps <- mean_total_var - mean_rep_var
  cat("Level:",level,"\n")
  cat("\tMean sample variance:",round(mean_total_var,3),"\n")
  cat("\tSE allotment approx",round(mean_total_var*0.9,3),"\n")   # 90%
  cat("\tPER allotment approx",round(mean_total_var*0.1,3),"\n")  # 10%
  #cat("\tNOISE allotment approx",round(mean_rep_var,3),"\n")
  #cat("\tSE allotment approx",round(mean_after_reps*0.9,3),"\n")   # 90%
  #cat("\tPER allotment approx",round(mean_after_reps*0.1,3),"\n")  # 10%
}









