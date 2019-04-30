library(RColorBrewer)
library(ggplot2)

source("include.R")

individuals <- c("DUI", "ECH", "LOG", "VET") # highly sampled
#individuals <- c("ACA", "DUI")
#individuals <- c("VEG") # 10-sample test individual

glom_data <- load_glommed_data(level="genus", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(filtered)

parsed_filename <- "original_data/parsed_enzymes.RData"
if(file.exists(parsed_filename)) {
  load("original_data/parsed_enzymes.RData")
} else {
  data.piphillin <- read_metagenomics(metadata, subset=FALSE)
  cat("Ordering Piphillin data...\n")
  data.piphillin <- intersect_order_metagenomics_samples(data.piphillin)
  cat("Saving Piphillin data...\n")
  save(data.piphillin, file="original_data/parsed_enzymes.RData")
}

metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)

for(indiv in individuals) {
  # plot matched 16S (genus-agglomerated) timecourse
  if(FALSE) {
    indiv_data <- subset_samples(glom_data, sname==indiv)
    family_data <- glom_counts(indiv_data, level="family", NArm=FALSE)
    filtered <- filter_data(family_data, count_threshold=3, sample_threshold=0.2)
    plot_timecourse_phyloseq(filtered, save_filename=paste(indiv, "_timecourse", sep=""), gapped=F, legend=F)
    plot_timecourse_phyloseq(filtered, save_filename=paste(indiv, "_timecourse_gapped", sep=""), gapped=T, legend=F)
  }
  
  # plot gapped metagenomics timecourse
  cat("Subsetting metagenomics data for",indiv,"\n")
  subset.piphillin <- subset_metagenomics_sname(data.piphillin, indiv, metadata.metagenomics)
  # this step takes a while (~1-2 min.) with the full enzyme list
  cat("Getting proportions from metagenomics data for",indiv,"\n")
  prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin, indiv, metadata.metagenomics)
  cat("Plotting metagenomics data for",indiv,"\n")
  plot_timecourse_metagenomics(prop.piphillin, save_filename=paste(indiv, "_metagenomics_timecourse", sep=""), gapped=F, legend=F)
  plot_timecourse_metagenomics(prop.piphillin, save_filename=paste(indiv, "_metagenomics_timecourse_gapped", sep=""), gapped=T, legend=F)
}
