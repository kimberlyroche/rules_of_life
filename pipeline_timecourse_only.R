source("include.R")

glom_data <- load_glommed_data(level="phylum", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66)
metadata <- read_metadata(filtered)

# ====================================================================================================================
# time courses for densely sampled individuals
# ====================================================================================================================

ACA_data <- subset_samples(filtered, sname=="ACA")
plot_timecourse_phyloseq(ACA_data, "ACA_phylum_timecourse", gapped=FALSE, legend=TRUE, legend_level="phylum")

PEB_data <- subset_samples(filtered, sname=="PEB")
plot_timecourse_phyloseq(PEB_data, "PEB_phylum_timecourse", gapped=FALSE, legend=TRUE, legend_level="phylum")
