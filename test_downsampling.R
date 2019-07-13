#library(stray)
devtools::load_all("/data/mukherjeelab/labraduck")

library(phyloseq)
library(dplyr)
library(driver)
library(tidyverse)

source("include.R")

# take the best-sampled individual and successively downsample them; do we see systematic effects
# of Sigma from increased uncertainty in the GP fits?

baboon <- "DUI"
level <- "family"
alr_ref <- 9

downsample_to <- seq(1, 0.2, by=-0.1)
downsample_to <- c(1, 0.2)
samples_per <- 100
downsample_reps <- 10

# read in and filter full data set at this phylogenetic level
glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)

# cut this down to the desired individual
indiv_data <- subset_samples(data, sname==baboon)

# get observations as differences from baseline in units of days
indiv_metadata <- read_metadata(indiv_data)
all_observations <- indiv_metadata$collection_date

Sigma_max <- NULL
Sigma_min <- NULL

for(ds in downsample_to) {

  Sigma <- array(NA, dim=c(26, 26, downsample_reps*samples_per))

  for(dsr in 1:downsample_reps) {

    Y <- otu_table(indiv_data)@.Data
    Y <- t(Y)

    # subsample
    ss <- ds*ncol(Y)
    ss_idx <- sort(sample(ncol(Y), ss))
    Y <- Y[,ss_idx]
    sub_observations <- all_observations[ss_idx]

    baseline_date <- sub_observations[1]
    observations <- sapply(sub_observations, function(x) round(difftime(x, baseline_date, units="days"))) + 1

    dim(observations) <- c(1, length(observations))
    colnames(Y) <- NULL
    rownames(Y) <- NULL

    D <- nrow(Y)
    N <- ncol(Y)

    cat("Y is",D,"x",N,"\n")

    dd_se <- 90
    dc <- 0.1 # desired SE minimum correlation
    se_sigma <- 1
    rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate

    period <- 365
    per_sigma <- 2
    rho_per <- 0.25

    se_weight <- 2
    per_weight <- 0.2

    Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
                         per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
                         (1e-7)*diag(ncol(X)) # pretty arbitrary
     
    # stray uses the D^th element as the ALR reference by default
    # do some row shuffling in Y to put the reference at the end
    if(!is.null(alr_ref)) {
      Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
    }
    
    # ALR prior covariance
    upsilon <- D-1+10 # lesser certainty
    GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference;
    # take diag as covariance over log abundances
    Xi <- GG%*%(diag(D))%*%t(GG)
    Xi <- Xi*(upsilon-D-1)
    
    alr_ys <- driver::alr((t(Y)+0.5), d=alr_ref)
    alr_means <- colMeans(alr_ys)
    Theta <- function(X) matrix(alr_means, D-1, ncol(X))
    
    fit <- stray::basset(Y, observations, upsilon, Theta, Gamma, Xi, n_samples=samples_per)
    fit.clr <- to_clr(fit)
    Sigma[,,(((dsr-1)*samples_per)+1):(dsr*samples_per)] <- fit.clr$Sigma[,,1:samples_per]
  }

  if(ds == max(downsample_to)) {
    Sigma_max <- Sigma
  }
  if(ds == min(downsample_to)) {
    Sigma_min <- Sigma
  }
  
  meanSigma <- apply(Sigma, c(1,2), mean)
  png(paste0("hist_",ds,".png"))
  plot(density(meanSigma))
  dev.off()
  png(paste0("meanSigma_",ds,".png"))
  image(meanSigma)
  dev.off()
  meanSigma <- cov2cor(meanSigma)
  png(paste0("hist_corr_",ds,".png"))
  plot(density(meanSigma))
  dev.off()
  png(paste0("meanSigma_corr_",ds,".png"))
  image(meanSigma)
  dev.off()

}

meanSmax <- apply(Sigma_max, c(1,2), mean)
meanSmin <- apply(Sigma_min, c(1,2), mean)
png(paste0("hist_maxmin.png"))
plot(density(meanSmax), xlim=c(min(c(meanSmax, meanSmin)), max(c(meanSmax, meanSmin))), ylim=c(0, 3.5))
lines(density(meanSmin), col="blue")
dev.off()
meanSmax <- cov2cor(meanSmax)
meanSmin <- cov2cor(meanSmin)
png(paste0("hist_maxmin_cor.png"))
plot(density(meanSmax), xlim=c(min(c(meanSmax, meanSmin)), max(c(meanSmax, meanSmin))), ylim=c(0, 2.5))
lines(density(meanSmin), col="blue")
dev.off()
