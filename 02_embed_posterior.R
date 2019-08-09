library(Rcpp)
library(ggplot2)
library(dplyr)
library(driver)
library(stray)
source("include.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Testing usage: Rscript 02_embed_posterior.R family Lambda", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
use_Riemann <- TRUE

# testing
#date_lower_limit <- "2001-10-01"
#date_upper_limit <- "2003-11-30"
date_lower_limit <- NULL
date_upper_limit <- NULL

sourceCpp("mat_dist.cpp")

indiv_obj <- fitted_individuals(level="family")
individuals <- indiv_obj$individuals
pattern_str <- indiv_obj$pattern_str
regexpr_str <- indiv_obj$regexpr_str

# get dimension D
fit_obj <- readRDS(paste0("subsetted_indiv_data/",level,"/",individuals[1],regexpr_str))
fit.clr <- to_clr(fit_obj$fit)
#Lambda <- fit_obj$fit$Lambda
Lambda <- fit.clr$Lambda
Sigma <- fit_obj$fit$Sigma

P <- dim(Lambda)[1]
N <- dim(Lambda)[2]
n_samples_subset <- 100

n_indiv <- length(individuals)
if(which_measure == "Sigma") {
  all_samples <- matrix(NA, P, P*n_samples_subset*n_indiv)
} else {
  # we'll use per-sample average to mitigate individuals having different N !
  all_samples <- matrix(NA, n_indiv*n_samples_subset, P)
}
indiv_labels <- c()
for(i in 1:n_indiv) {
  fn <- paste0("subsetted_indiv_data/",level,"/",individuals[i],regexpr_str)
  fit <- readRDS(fn)$fit
  # to ILR
  #V <- driver::create_default_ilr_base(ncategories(fit))
  #fit.ilr <- to_ilr(fit, V)
  #Lambda <- fit.ilr$Lambda
  #Sigma <- fit.ilr$Sigma
  # to CLR
  fit.clr <- to_clr(fit)
  Lambda <- fit.clr$Lambda
  Sigma <- fit.clr$Sigma

  if(which_measure == "Sigma") {
    Sigma <- Sigma[,,1:n_samples_subset]
    all_samples[,((i-1)*P*n_samples_subset+1):(i*P*n_samples_subset)] <- Sigma
    indiv_labels <- c(indiv_labels, rep(individuals[i], n_samples_subset))
  } else {
    collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) })) # 100 x P
    all_samples[((i-1)*n_samples_subset+1):(i*n_samples_subset),] <- collLambda
    indiv_labels <- c(indiv_labels, rep(individuals[i], n_samples_subset))
  }
}

if(which_measure == "Sigma") {
  distance_mat <- mat_dist(all_samples, n_indiv, n_samples_subset)
} else {
  distance_mat <- dist(all_samples)
}

fit <- cmdscale(distance_mat, eig=TRUE, k=3) # k is the number of dim
cat("Eigval 1:",fit$eig[1],"\n")
cat("Eigval 2:",fit$eig[2],"\n")
cat("Eigval 3:",fit$eig[3],"\n")
cat("Eigval 4:",fit$eig[4],"\n")

df <- data.frame(x=fit$points[,1], y=fit$points[,2], z=fit$points[,3], labels=indiv_labels)
saveRDS(df, paste0("plots/basset/",level,"/",which_measure,"_ordination.rds"))

# centroids are useful for labeling plots
df_centroids <- df %>% group_by(labels) %>% summarise(mean_x=mean(x), mean_y=mean(y), mean_z=mean(z))
saveRDS(df_centroids, paste0("plots/basset/",level,"/",which_measure,"_ordination_centroids.rds"))

