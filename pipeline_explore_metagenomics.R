library(RColorBrewer)
library(ggplot2)

source("include.R")

data.16S <- read_data()
md.16S <- read_metadata(data.16S)

data.piphillin <- read_metagenomics(md.16S)
subset.piphillin <- subset_metagenomics_sname(data.piphillin, "ACA", md.16S)
prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin, "ACA", md.16S)
plot_metagenomics_timecourse(prop.piphillin)
