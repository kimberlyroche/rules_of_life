source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for task missing!\n")
}

task <- as.numeric(args[1])

if(task == 1) {
  # agglomerate taxa
  perform_agglomeration(level="species") # original data
  #perform_agglomeration(level="genus", replicates=TRUE) # data with replicates
  #perform_agglomeration(level="species", replicates=TRUE) # data with replicates
}

if(task == 2) {
  # render histograms of sample frequency, etc.
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  histogram_indiv_samples(filtered)
  histogram_sample_density(filtered, "weeks")
}

if(task == 3) {
  # plot sample time courses
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  perform_mult_timecourse(filtered, c("ACA", "DUI", "CAI", "COB", "DAS"))
}

if(task == 4) {
  # plot weekly autocorrelation out to 26 weeks; use zero-filter stringency of 20%
  t <- 0.2
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=t)
  lags <- calc_autocorrelation(filtered, lag.max=26, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk_",round(100*t),sep=""))
}

if(task == 5) {
  # plot weekly autocorrelation out to 26 weeks; use zero-filter stringency of 66%
  t <- 0.66
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=t)
  lags <- calc_autocorrelation(filtered, lag.max=26, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk_",round(100*t),sep=""))
}

if(task == 6) {
  # plot weekly autocorrelation out to 26 weeks; use zero-filter stringency of 90%
  t <- 0.9
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=t)
  lags <- calc_autocorrelation(filtered, lag.max=26, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk_",round(100*t),sep=""))
}

if(task == 7) {
  # plot weekly autocorrelation out to 52 weeks
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  lags <- calc_autocorrelation(filtered, lag.max=52, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_52wk")
}

if(task == 8) {
  # plot monthly autocorrelation out to 12 months
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  lags <- calc_autocorrelation(filtered, lag.max=12, resample=TRUE, date_diff_units="months")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_12mo")
}

if(task == 9) {
  # plot monthly autocorrelation out to 24 months
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  lags <- calc_autocorrelation(filtered, lag.max=24, resample=TRUE, date_diff_units="months", resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_24mo")
}

if(task == 10) {
  # plot monthly autocorrelation out to 36 months
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  lags <- calc_autocorrelation(filtered, lag.max=36, resample=TRUE, date_diff_units="months", resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_36mo")
}

if(task == 11) {
  # plot seasonal autocorrelation out to 5 years (~lag.max == 1)
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  lags <- calc_autocorrelation(filtered, lag.max=11, resample=TRUE, date_diff_units="seasons")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_11season")
}

if(task == 12) {
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
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

if(task == 13) {
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

if(task == 14) {
  # get distribution of within or between group correlation
  load("glom_data_species.RData")
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)

  log_ratios <- apply_ilr(filtered)
  log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))

  md <- read_metadata(filtered)
  groups <- unique(md$grp)
  within_corr <- c()
  between_corr <- c()
  for(i in 1:(length(groups)-1)) {
    # quick and dirty: make sure group has > some minimum number of samples
    if(nsamples(subset_samples(filtered, grp==groups[i])) > 100) {
      for(j in i:length(groups)) {
        if(nsamples(subset_samples(filtered, grp==groups[j])) > 100) {
          for(k in 1:100) {
            # get sample 1
            sample_set.1 <- subset_samples(filtered, grp==groups[i])
            idx <- sample(nsamples(sample_set.1))[1]
            sample_idx.1 <- read_metadata(sample_set.1)$sample_id[idx]
            sample.1 <- subset_samples(sample_set.1, sample_id==sample_idx.1)
            sample_season.1 <- read_metadata(sample.1)$season
            sample_sname.1 <- read_metadata(sample.1)$sname

            # get sample 2
            sample_set.2 <- subset_samples(filtered, grp==groups[j])
            # match season
            sample_set.2 <- subset_samples(sample_set.2, season == sample_season.1)
            # exclude same-individual
            sample_set.2 <- subset_samples(sample_set.2, sname != sample_sname.1)
            idx <- sample(nsamples(sample_set.2))[1]
            sample_idx.2 <- read_metadata(sample_set.2)$sample_id[idx]
            sample.2 <- subset_samples(sample_set.2, sample_id==sample_idx.2)

            sample_lr.1 <- log_ratios[,colnames(log_ratios) %in% sample_idx.1]
            sample_lr.2 <- log_ratios[,colnames(log_ratios) %in% sample_idx.2]

            y.t <- as.vector(sample_lr.1)
            y.tt <- sqrt(y.t%*%y.t)
            y.h <- as.vector(sample_lr.2)
            y.hh <- sqrt(y.h%*%y.h)
            corr_est <- (y.t%*%y.h)/(y.tt*y.hh)
            if(i == j) {
              within_corr[length(within_corr)+1] <- corr_est
              cat("Within",groups[i],":",corr_est,"\n")
            } else {
              between_corr[length(between_corr)+1] <- corr_est
              cat("Between",groups[i],"and",groups[j],":",corr_est,"\n")
            }
          }
        }
      }
    }
  }
  cat("Average within:",mean(within_corr),"\n")
  cat("Average between:",mean(between_corr),"\n")
}
