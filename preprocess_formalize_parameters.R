# get empirical estimates for variance to use in GP kernels

source("include/general.R")

# PHYLUM level estimates

levels <- c("phylum", "family")

for(level in levels) {
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, collapse_level=level)

  metadata <- sample_data(data)
  indiv_sd <- c()
  for(indiv in unique(metadata$sname)) {
    cat("Parsing individual",indiv,"\n")
    data_subset <- subset_samples(data, sname == indiv)
    if(nsamples(data_subset) > 1) {
      alr_data <- driver::alr(otu_table(data_subset)@.Data + 0.5)
      tax_sd <- apply(alr_data, 2, sd)
      indiv_sd <- c(indiv_sd, tax_sd)
    }
  }

  # get average ALR travel: min-max difference within a coordinate within an individual
  mean_sd <- mean(indiv_sd)
  cat("Level:",level,"\n")
  cat("\tMean ALR deviation:",round(mean_sd,3),"\n")
  cat("\tSE allotment approx",round(mean_sd*0.9,3),"\n")   # 90%
  cat("\tPER allotment approx",round(mean_sd*0.1,3),"\n")  # 10%
}
