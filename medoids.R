library(cluster)
source("include.R")

filtered <- filter_data()

sample_data <- subset_samples(filtered, sname %in% c("DUI", "ACA"))
sample_ilr <- apply_ilr(sample_data)
sample_ilr <- t(apply(samples_ilr, 1, function(x) x - mean(x)))

diss_mat <- 1-cov2cor(cov(sample_ilr))

test <- pam(diss_mat, 10, metric="euclidean", do.swap = TRUE)
