source("include.R")

#args = commandArgs(trailingOnly=TRUE)
#if(length(args)==0) {
#  stop("Argument for task missing!\n")
#}

#task <- as.numeric(args[1])
task <- 8

if(task == 1) {
  # agglomerate to genus
  perform_agglomeration()
}

#if(!exists("filtered")) {
#  filtered <- filter_data()
#  #filtered_subset <- subset_samples(filtered, sname %in% best_sampled)
#  #filtered_subset <- subset_samples(filtered, sname %in% over_50)
#}

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
  # plot weekly autocorrelation out to 26 weeks
  lags <- calc_autocorrelation(filtered, lag.max=26, resample=TRUE, resample_rate=0.5)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_26wk")
}

if(task == 5) {
  # plot weekly autocorrelation out to 52 weeks
  lags <- calc_autocorrelation(filtered, lag.max=52, resample=TRUE, resample_rate=0.66)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_52wk")
}

if(task == 6) {
  # plot monthly autocorrelation out to 12 months
  lags <- calc_autocorrelation(filtered, lag.max=12, resample=TRUE, date_diff_units="months")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_12mo")
}

if(task == 7) {
  # plot monthly autocorrelation out to 24 months
  lags <- calc_autocorrelation(filtered, lag.max=24, resample=TRUE, date_diff_units="months", resample_rate=0.66)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_24mo")
}

if(task == 8) {
  # plot monthly autocorrelation out to 36 months
  lags <- calc_autocorrelation(filtered, lag.max=36, resample=TRUE, date_diff_units="months", resample_rate=0.66)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_36mo")
}

if(task == 9) {
  # plot scrambled correlation matrices using grouped sample
  # batch stuff
  visualize_groupwise_covariance(filtered, "plate", sample=50)
  visualize_groupwise_covariance(filtered, "flow_cell_lane", sample=1000)
  visualize_groupwise_covariance(filtered, "library_pool", sample=200)
  visualize_groupwise_covariance(filtered, "extract_dna_conc_ng", sample=100)
  # biological stuff
  visualize_groupwise_covariance(filtered, "sex", sample=1000)
  visualize_groupwise_covariance(filtered, "season", sample=1000)
  visualize_groupwise_covariance(filtered, "grp", sample=100)
  visualize_groupwise_covariance(filtered, "matgrp", sample=100)
  visualize_groupwise_covariance(filtered, "sname", sample=100)
  visualize_groupwise_covariance(filtered, "age", sample=500)
}

if(task == 10) {
  # subset to a manageable sample size
  md <- read_metadata(filtered)
  sample_ids <- md$sample_id[sample(nsamples(filtered))[1:5000]]
  baboon_counts <- subset_samples(filtered, sample_id %in% sample_ids)
  # fix these, so they don't plot as continuous values
  sample_data(baboon_counts)$grp <- as.factor(sample_data(baboon_counts)$grp)
  sample_data(baboon_counts)$matgrp <- as.factor(sample_data(baboon_counts)$matgrp)
  sample_data(baboon_counts)$age <- cut(sample_data(baboon_counts)$age,
                                        breaks=c(0, 4.5, 19, Inf), labels=c("juvenile", "adult", "geriatric"))

  # PCoA ordination using Bray-Curtis distance
  ord <- ordinate(baboon_counts, "PCoA", "bray")

  cat("Plotting seasonal labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="season")
  ggsave("plots/season_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting collection group-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="grp")
  ggsave("plots/grp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting sex-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="sex")
  ggsave("plots/sex_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting maternal group-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="matgrp")
  ggsave("plots/matgrp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
  cat("Plotting maternal age-wise labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="age")
  ggsave("plots/age_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

  baboon_counts <- subset_samples(filtered, sname %in% best_sampled)

  # PCoA ordination using Bray-Curtis distance
  ord <- ordinate(baboon_counts, "PCoA", "bray")

  cat("Plotting individual labeling on ordination...\n")
  p <- plot_ordination(baboon_counts, ord, type="sample", color="sname")
  ggsave("plots/sname_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
}
