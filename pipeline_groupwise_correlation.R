# estimates within-group and between-group correlation

source("include.R")

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)

# apply log ratio transform
log_ratios <- apply_ilr(filtered)
log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))

md <- read_metadata(filtered)
groups <- unique(md$grp)
within_corr <- c()
between_corr <- c()
for(i in 1:(length(groups)-1)) {
  # quick and dirty: skip groups for pairwise comparison if each doesn't have
  # some minimum number of samples (e.g. 100)
  if(nsamples(subset_samples(filtered, grp==groups[i])) > 100) {
    for(j in i:length(groups)) {
      if(nsamples(subset_samples(filtered, grp==groups[j])) > 100) {
        for(k in 1:100) {
          # sample group A
          sample_set.1 <- subset_samples(filtered, grp==groups[i])
          idx <- sample(nsamples(sample_set.1))[1]
          sample_idx.1 <- read_metadata(sample_set.1)$sample_id[idx]
          sample.1 <- subset_samples(sample_set.1, sample_id==sample_idx.1)
          sample_season.1 <- read_metadata(sample.1)$season
          sample_sname.1 <- read_metadata(sample.1)$sname

          # sample group B
          sample_set.2 <- subset_samples(filtered, grp==groups[j])
          # match season
          sample_set.2 <- subset_samples(sample_set.2, season == sample_season.1)
          # exclude samples from the same individual
          sample_set.2 <- subset_samples(sample_set.2, sname != sample_sname.1)
          idx <- sample(nsamples(sample_set.2))[1]
          sample_idx.2 <- read_metadata(sample_set.2)$sample_id[idx]
          sample.2 <- subset_samples(sample_set.2, sample_id==sample_idx.2)

          sample_lr.1 <- log_ratios[,colnames(log_ratios) %in% sample_idx.1]
          sample_lr.2 <- log_ratios[,colnames(log_ratios) %in% sample_idx.2]

          # calculate "correlation" (cosine angle)
          y.t <- as.vector(sample_lr.1)
          y.tt <- sqrt(y.t%*%y.t)
          y.h <- as.vector(sample_lr.2)
          y.hh <- sqrt(y.h%*%y.h)
          corr_est <- (y.t%*%y.h)/(y.tt*y.hh)
          if(i == j) {
            within_corr[length(within_corr)+1] <- corr_est
          } else {
            between_corr[length(between_corr)+1] <- corr_est
          }
        }
      }
    }
  }
}

cat("Average within-group correlation:",mean(within_corr),"\n")
cat("Average between-group correlation:",mean(between_corr),"\n")
