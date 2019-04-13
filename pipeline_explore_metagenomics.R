library(RColorBrewer)
library(ggplot2)

source("include.R")

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
md.16S <- read_metadata(filtered)

indiv <- "COO"

# plot gapped metagenomics timecourse
data.piphillin <- read_metagenomics(md.16S)
subset.piphillin <- subset_metagenomics_sname(data.piphillin, indiv, md.16S)
prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin, indiv, md.16S)
plot_timecourse_metagenomics(prop.piphillin, save_filename=paste(indiv, "_metagenomics_timecourse", sep=""), gapped=T, legend=F)

# plot matched 16S (genus-agglomerated) timecourse
perform_mult_timecourse(filtered, c(indiv), gapped=T, legend=F)
