source("include/R/general.R")
source("include/R/metagenomics")

individuals <- c("DUI", "ECH", "LOG", "VET") # highly sampled

level <- "genus"
data <- load_and_filter(level)
metadata <- read_metadata(data)

parsed_filename <- paste0(data_dir,"Piphillin_20190222/parsed_enzymes.RData")
if(file.exists(parsed_filename)) {
  load(parsed_filename)
} else {
  data.piphillin <- read_metagenomics(metadata, subset=FALSE)
  cat("Ordering Piphillin data...\n")
  data.piphillin <- intersect_order_metagenomics_samples(data.piphillin)
  cat("Saving Piphillin data...\n")
  save(data.piphillin, file=paste0(data_dir,"parsed_enzymes.RData")
}

metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)

for(indiv in individuals) {
  # plot matched 16S (genus-agglomerated) timecourse
  if(FALSE) {
    indiv_data <- subset_samples(glom_data, sname==indiv)
    family_data <- glom_counts(indiv_data, level="family", NArm=FALSE)
    filtered <- filter_data(family_data, level=level)
    plot_timecourse_phyloseq(filtered, save_filename=paste0(plot_dir,indiv,"_timecourse"), gapped=F, legend=F)
    plot_timecourse_phyloseq(filtered, save_filename=paste0(plot_dir,indiv, "_timecourse_gapped"), gapped=T, legend=F)
  }
  
  # plot gapped metagenomics timecourse
  cat("Subsetting metagenomics data for",indiv,"\n")
  subset.piphillin <- subset_metagenomics_sname(data.piphillin, indiv, metadata.metagenomics)
  
  # find the large-proportion enzymes
  subset.prop <- prop.table(subset.piphillin, 2)
  subset.mean <- rowMeans(subset.prop)
  temp <- subset.mean[order(subset.mean)]
  cat("Top 10 (proportional) ",indiv,":\n",sep="")
  cat(names(temp[(length(temp)-10):length(temp)]),sep="\n")
  
  if(FALSE) {
    # this step takes a while (~5 min.) with the full enzyme list
    cat("Getting proportions from metagenomics data for",indiv,"\n")
    prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin, indiv, metadata.metagenomics)
    cat("Plotting metagenomics timecourse data for",indiv,"\n")
    plot_timecourse_metagenomics(prop.piphillin, save_filename=paste0(plot_dir,indiv,"_metagenomics_timecourse"), gapped=F, legend=F)
    plot_timecourse_metagenomics(prop.piphillin, save_filename=paste0(plot_dir,indiv,"_metagenomics_timecourse_gapped"), gapped=T, legend=F)
  }
}
