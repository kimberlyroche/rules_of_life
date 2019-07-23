source("include.R")

glom_data <- load_glommed_data(level="phylum", replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)

# get average ALR "travel" (min-max delta within a coordinate within an individual)
cat("Parsing individuals...\n")
metadata <- sample_data(data)
indiv_sd <- c()
for(indiv in unique(metadata$sname)) {
  data_subset <- subset_samples(data, sname == indiv)
  if(nsamples(data_subset) > 1) {
    alr_data <- driver::alr(otu_table(data_subset)@.Data + 0.5)
    tax_sd <- apply(alr_data, 2, sd)
    indiv_sd <- c(indiv_sd, tax_sd)
  }
}
mean_sd <- mean(indiv_sd)
cat("PHYLUM\n")
cat("\tMean ALR deviation:",mean_sd,"\n")
cat("\tSE allotment approx",mean_sd*0.9,"\n")
cat("\tPER allotment approx",mean_sd*0.1,"\n")

glom_data <- load_glommed_data(level="family", replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)

# get average ALR "travel" (min-max delta within a coordinate within an individual)
cat("Parsing individuals...\n")
metadata <- sample_data(data)
indiv_sd <- c()
for(indiv in unique(metadata$sname)) {
  data_subset <- subset_samples(data, sname == indiv)
  if(nsamples(data_subset) > 1) {
    alr_data <- driver::alr(otu_table(data_subset)@.Data + 0.5)
    tax_sd <- apply(alr_data, 2, sd)
    indiv_sd <- c(indiv_sd, tax_sd)
  }
}
mean_sd <- mean(indiv_sd)
cat("FAMILY\n")
cat("\tMean ALR deviation:",mean_sd,"\n")
cat("\tSE allotment approx",mean_sd*0.9,"\n")
cat("\tPER allotment approx",mean_sd*0.1,"\n")

