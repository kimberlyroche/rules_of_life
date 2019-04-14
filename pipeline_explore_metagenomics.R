library(RColorBrewer)
library(ggplot2)

source("include.R")

individuals <- c("ACA", "DUI")

individuals <- c("VEG") # 10-sample test individual

for(indiv in individuals) {
  
  # plot matched 16S (genus-agglomerated) timecourse
  
  glom_data <- load_glommed_data(level="genus", replicates=TRUE)
  indiv_data <- subset_samples(glom_data, sname==indiv)
  genus_data <- glom_counts(indiv_data, level="family", NArm=FALSE)
  filtered <- filter_data(genus_data, count_threshold=3, sample_threshold=0.2)
  plot_timecourse_phyloseq(filtered, save_filename=paste(indiv, "_timecourse", sep=""), gapped=F, legend=F)
  plot_timecourse_phyloseq(filtered, save_filename=paste(indiv, "_timecourse_gapped", sep=""), gapped=T, legend=F)
  
  # plot gapped metagenomics timecourse
  
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  metadata <- read_metadata(filtered)
  
  data.piphillin <- read_metagenomics(metadata)
  metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)
  subset.piphillin <- subset_metagenomics_sname(data.piphillin, indiv, metadata.metagenomics)
  prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin, indiv, metadata.metagenomics)
  plot_timecourse_metagenomics(prop.piphillin, save_filename=paste(indiv, "_metagenomics_timecourse", sep=""), gapped=F, legend=F)
  plot_timecourse_metagenomics(prop.piphillin, save_filename=paste(indiv, "_metagenomics_timecourse_gapped", sep=""), gapped=T, legend=F)

  }