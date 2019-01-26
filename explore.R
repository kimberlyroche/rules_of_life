source("include.R")

task <- 6

if(task == 1) {
  # agglomerate to genus
  perform_agglomeration()
}

if(!exists("filtered")) {
  filtered <- filter_data()
}

if(task == 2) {
  # render histograms of sample frequency, etc.
  histogram_indiv_samples(filtered)
  histogram_sample_density(filtered, "weeks")
}

if(task == 3) {
  # plot sample time courses
  perform_mult_timecourse(filtered, c("ACA", "DUI", "CAI", "COB", "DAS"))
}

if(task == 4) {
  # plot (weekly) autocorrelation
  plot_autocorrelation(filtered)
}

if(task == 5) {
  # plot scrambled correlation matrices using grouped sample
  visualize_groupwise_covariance(filtered, "sex", sample=1000)
  visualize_groupwise_covariance(filtered, "season", sample=1000)
  visualize_groupwise_covariance(filtered, "grp", sample=100)
  visualize_groupwise_covariance(filtered, "matgrp", sample=100)
  visualize_groupwise_covariance(filtered, "sname", sample=100)
  visualize_groupwise_covariance(filtered, "age", sample=500)
}

if(task == 6) {
  # subset to a manageable sample size
  md <- read_metadata(filtered)
  sample_ids <- md$sample_id[sample(nsamples(filtered))[1:2000]]
  baboon_counts <- subset_samples(filtered, sample_id %in% sample_ids)
  # fix these, so they don't plot as continuous values
  sample_data(baboon_counts)$grp <- as.factor(sample_data(baboon_counts)$grp)
  sample_data(baboon_counts)$matgrp <- as.factor(sample_data(baboon_counts)$matgrp)

  # PCoA ordination using Bray-Curtis distance
  ord <- ordinate(baboon_counts, "PCoA", "bray")

  cat("Plotting seasonal labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="season")
  ggsave("season_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting collection group-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="grp")
  ggsave("grp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting sex-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="sex")
  ggsave("sex_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting maternal group-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="matgrp")
  ggsave("matgrp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

  baboon_counts <- subset_samples(filtered, sname %in% best_sampled)

  # PCoA ordination using Bray-Curtis distance
  ord <- ordinate(baboon_counts, "PCoA", "bray")

  cat("Plotting individual labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="sname")
  ggsave("sname_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
}
