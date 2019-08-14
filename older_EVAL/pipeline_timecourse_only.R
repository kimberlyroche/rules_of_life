source("include.R")

glom_data <- load_glommed_data(level="phylum", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66)
metadata <- read_metadata(filtered)

# ====================================================================================================================
# time courses for densely sampled individuals
# ====================================================================================================================

individuals <- c("DUN", "ZIZ")
for(indiv in individuals) {
  indiv_data <- subset_samples(filtered, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(indiv,"_phylum_timecourse"), gapped=FALSE, legend=TRUE, legend_level="phylum")
}
