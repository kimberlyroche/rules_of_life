source("include/R/general.R")
source("include/R/GP.R")

levels <- c("phylum", "family", "genus")
levels <- c("family", "genus")
count_threshold <- 1
sample_threshold <- 0.05

# -------------------------------------------------------------------------------------------------------
# WHAT PROPORTION OF TAXA GET THROWN INTO NA FOR A HIGH ALPHA-DIV INDIVIDUAL?
# -------------------------------------------------------------------------------------------------------

cat("Testing presence between hosts...\n")
for(level in levels) {
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  glom_data_subset <- subset_samples(glom_data, sname %in% over_50)
  data <- filter_data(glom_data_subset, count_threshold=count_threshold, sample_threshold=sample_threshold, collapse_level=level, verbose=TRUE)
  cat("^^^Assuming index of OTHER is 2. CHECK THIS.\n")
  for(indiv in c("ZIZ", "GAN")) {
    indiv_data <- subset_samples(data, sname==indiv)
    counts <- otu_table(indiv_data)@.Data # taxa are columns
    prop_NA <- apply(counts, 1, function(x) x[2]/sum(x))
    names(prop_NA) <- NULL
    cat("Mean, sd proportion of NA in",indiv,":", round(mean(prop_NA), 3),",", round(sd(prop_NA), 3),"\n")
  }
}

# -------------------------------------------------------------------------------------------------------
# WHAT PROPORTION OF TAXA ARE PRESENT ACROSS HOSTS?
# -------------------------------------------------------------------------------------------------------

if(FALSE) {
cat("Testing presence between hosts...\n")
for(level in levels) {
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, count_threshold=count_threshold, sample_threshold=sample_threshold, collapse_level=level, verbose=TRUE)
  # what proportion of taxa are present (non-zero) in at least 10% of all individual's samples?
  taxa <- numeric(ntaxa(data))
  for(i in 1:length(over_50)) {
    data2 <- subset_samples(data, sname==over_50[i])
    threshold <- 1
    counts <- otu_table(data2)@.Data
    taxa <- taxa + as.numeric(as.vector(unlist(apply(counts, 2, function(x) { sum(x != 0) > threshold }))))
  }
  cat("Level:",level,"\n")
  print(table(taxa))
}
}

# -------------------------------------------------------------------------------------------------------
# WHAT PROPORTION OF TAXA ARE PRESENT IN >= 10% SAMPLES  ACROSS HOSTS?
# -------------------------------------------------------------------------------------------------------

if(FALSE) {
cat("Testing >10% presence between hosts...\n")
for(level in levels) {
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.2, collapse_level=level, verbose=FALSE)
  # what proportion of taxa are present (non-zero) in at least 10% of all individual's samples?
  taxa <- numeric(ntaxa(data))
  for(i in 1:length(over_50)) {
    data2 <- subset_samples(data, sname==over_50[i])
    threshold <- 0.1*length(phyloseq::nsamples(data2))
    counts <- otu_table(data2)@.Data
    taxa <- taxa + as.numeric(as.vector(unlist(apply(counts, 2, function(x) { sum(x != 0) > threshold }))))
  }
  cat("Level:",level,"\n")
  print(table(taxa))
}
}

# -------------------------------------------------------------------------------------------------------
# DOES DISTANCE OVER BASELINE CORRELATE WITH DISTANCE OVER COVARIANCE MATRICES?
# -------------------------------------------------------------------------------------------------------

if(FALSE) {
cat("Testing >10% presence between hosts...\n")
for(level in levels) {
  cat("Level:",level,"\n")

  if(level == "phylum") {
    se_weight <- sqrt(1.956)
    per_weight <- sqrt(0.217)
    wn_weight <- sqrt(0.233)
  }
  if(level == "family") {
    se_weight <- sqrt(2.090)
    per_weight <- sqrt(0.232)
    wn_weight <- sqrt(0.381)
  }
  if(level == "genus") {
    se_weight <- sqrt(2.048)
    per_weight <- sqrt(0.228)
    wn_weight <- sqrt(0.416)
  }

  baboons <- over_50
  all_samples_Sigma <- NULL
  all_samples_Lambda <- NULL
  P <- NULL
  V <- NULL
  for(i in 1:length(baboons)) {
    baboon <- baboons[i]
    fit_obj <- fit_GP(baboon, level, se_weight=se_weight, per_weight=per_weight, wn_weight=wn_weight,
                      dd_se=90, save_append="", date_lower_limit=NULL, date_upper_limit=NULL, verbose=FALSE, mean_only=TRUE)
    dim(fit_obj$fit$Eta) <- c(nrow(fit_obj$fit$Eta), ncol(fit_obj$fit$Eta), 1)
    dim(fit_obj$fit$Lambda) <- c(nrow(fit_obj$fit$Lambda), ncol(fit_obj$fit$Lambda), 1)
    dim(fit_obj$fit$Sigma) <- c(nrow(fit_obj$fit$Sigma), ncol(fit_obj$fit$Sigma), 1)
    if(is.null(V)) {
      V <- driver::create_default_ilr_base(ncategories(fit_obj$fit))
    }
    fit.ilr <- to_ilr(fit_obj$fit, V)
    Lambda <- fit.ilr$Lambda
    Sigma <- fit.ilr$Sigma
     if(is.null(all_samples_Sigma)) {
      P <- fit.ilr$D-1
    }
  
    # note: this is calculated incorrectly in PibbleCollapse_Uncollapse.cpp
    upsilonN <- fit.ilr$upsilon + fit.ilr$N
    Sigma <- Sigma / (upsilonN - fit.ilr$D)
    Sigma <- Sigma / (upsilonN - P - 1)
    # this should be the correct invWishart mean
    if(is.null(all_samples_Sigma)) {
      all_samples_Sigma <- matrix(NA, P, P*length(baboons))
    }
    all_samples_Sigma[,((i-1)*P+1):(i*P)] <- Sigma
    if(is.null(all_samples_Lambda)) {
      all_samples_Lambda <- matrix(NA, length(baboons), P)
    }
    collLambda <- apply(Lambda[,,1], 1, mean) # P x 1
    all_samples_Lambda[i,] <- collLambda
  }
  distance_mat <- Riemann_dist_samples(all_samples_Sigma, length(baboons), 1)
  distances_Sigma <- distance_mat[lower.tri(distance_mat)]
  distances_Lambda <- c(dist(all_samples_Lambda))

  cat("Correlation between distances:",cor(distances_Sigma, distances_Lambda),"\n")

  df <- data.frame(x=distances_Lambda, y=distances_Sigma)
  p <- ggplot(df) +
       geom_point(aes(x=x, y=y)) +
       xlab("distance(baseline)") +
       ylab("distance(covariance)")
  ggsave(paste0(plot_dir,"posterior_dist_vs_dist_",level,".png"), plot=p, scale=1.5, width=5, height=5, units="in", dpi=100)
}
}
