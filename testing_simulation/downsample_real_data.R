# this file repeatedly downsamples real data from the best sampled individual and fits the model,
# helping us visualize the effect of sample number on the quality of posterior estimates

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

sourceCpp(file.path(relative_path,"include/Rcpp/Riemann_dist.cpp"))

# calculate distances over a bunch of Sigma samples and return as vector
# Sigma_samples is a D x D x N array
get_Sigma_spread <- function(Sigma_samples) {
  d <- c()
  dcor <- c()
  N <- dim(Sigma_samples)[3]
  # indexing is because we're calculating lower triangular only
  for(i in 1:(N-1)) {
    # repeat for covariance estimate and correlation estimate of Sigma
    A <- Sigma_samples[,,i]
    Acor <- cov2cor(A)
    for(j in (i+1):N) {
      B <- Sigma_samples[,,j]
      Bcor <- cov2cor(B)
      d <- c(d, Riemann_dist_pair(A, B))
      dcor <- c(dcor, Riemann_dist_pair(Acor, Bcor))
    }
  }
  return(list(distances=d, distances_corr=dcor))
}

baboon <- "DUI"
level <- "family"
alr_ref <- NULL

downsample_to <- seq(1, 0.2, by=-0.1)
samples_per <- 100 # posterior sample number
downsample_reps <- 10

data <- load_and_filter(level)
indiv_data <- subset_samples(data, sname==baboon)
indiv_metadata <- read_metadata(indiv_data)
all_observations <- indiv_metadata$collection_date
D <- ntaxa(indiv_data)

# these are whole Sigma sample sets associated with
#   max sample percent (100%)
#   median sample percent
#   min sample percent (20%)
# we'll use these to crudely ask about how posterior spread is changing
Sigma_max <- NULL
Sigma_med <- NULL
Sigma_min <- NULL

# iterate down from 100% of samples to 20% of samples
# for no particular reason, we're just going to stick to the ALR
for(ds in downsample_to) {
  coord_dim <- D # ALR
  # Sigma samples concatenated horizontally
  Sigma <- array(NA, dim=c(D-1, D-1, downsample_reps*samples_per))

  # repeat downsampling a few times
  for(dsr in 1:downsample_reps) {
    # Y is taxa (rows) x samples (columns)
    Y <- t(otu_table(indiv_data)@.Data)
    # subsample the counts
    ss <- ds*ncol(Y)
    ss_idx <- sort(sample(ncol(Y), ss))
    Y <- Y[,ss_idx]
    sub_observations <- all_observations[ss_idx]
    # subsample the observations to match and transform them into 1-indexed timepoints (days)
    baseline_date <- sub_observations[1]
    observations <- sapply(sub_observations, function(x) round(difftime(x, baseline_date, units="days"))) + 1
    dim(observations) <- c(1, length(observations))
    colnames(Y) <- NULL; rownames(Y) <- NULL

    N <- ncol(Y)

    # use a vanilla parameterization for basset
    dd_se <- 90
    dc <- 0.1
    se_sigma <- 1
    rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate
    period <- 365
    per_sigma <- 2
    rho_per <- 0.25
    se_weight <- 2
    per_weight <- 0.2

    Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
                         per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
                         (1e-8)*diag(ncol(X)) # pretty arbitrary
     
    if(!is.null(alr_ref)) {
      Y <- reorient_count_matrix(Y, alr_ref)
    }
    
    prior_obj <- default_ALR_prior(D)    
    alr_ys <- driver::alr((t(Y)+pc), d=alr_ref)
    alr_means <- colMeans(alr_ys)
    Theta <- function(X) matrix(alr_means, D-1, ncol(X))
    
    fit <- stray::basset(Y, observations, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=samples_per)
    Sigma[,,(((dsr-1)*samples_per)+1):(dsr*samples_per)] <- fit$Sigma[,,1:samples_per]
  }

  if(ds == max(downsample_to)) {
    Sigma_max <- Sigma
  }
  if(ds == median(downsample_to)) {
    Sigma_med <- Sigma
  }
  if(ds == min(downsample_to)) {
    Sigma_min <- Sigma
  }

  # take the element-wise mean over all these samples and visualize it  
  meanSigma <- apply(Sigma, c(1,2), mean)
  png(file.path(relative_path,plot_dir,paste0("downsampled_meanSigma_",ds,".png")))
  image(meanSigma)
  dev.off()
  # repeat as correlation
  meanSigma <- cov2cor(meanSigma)
  png(file.path(relative_path,plot_dir,paste0("downsampled_meanSigma_corr_",ds,".png")))
  image(meanSigma)
  dev.off()
}

# plot Riemannian distance for MAX, MEDIAN, and MIN sample density (crude)
d_max <- get_Sigma_spread(Sigma_max)
d_med <- get_Sigma_spread(Sigma_med)
d_min <- get_Sigma_spread(Sigma_min)
all_d <- c(d_max$distances, d_med$distances, d_min$distances)
xlim <- c(min(all_d), max(all_d))
png(file.path(relative_path,plot_dir,paste0("downsampled_distance_distribution.png")))
plot(density(d_max$distances), xlim=xlim)
lines(density(d_med$distances), col="red")
lines(density(d_min$distances), col="blue")
dev.off()

# repeat as correlation
all_d <- c(d_max$distances_corr, d_med$distances_corr, d_min$distances_corr)
xlim <- c(min(all_d), max(all_d))
png(file.path(relative_path,plot_dir,"downsampled_distance_distribution_corr.png"))
plot(density(d_max$distances_corr), xlim=xlim)
lines(density(d_med$distances_corr), col="red")
lines(density(d_min$distances_corr), col="blue")
dev.off()
