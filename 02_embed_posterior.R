library(Rcpp)
library(ggplot2)
library(dplyr)
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

sourceCpp("cov_viz_test.cpp")

pattern_str <- "*_bassetfit.RData"
regexpr_str <- "_bassetfit.RData"
fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern=pattern_str, full.names=TRUE, recursive=FALSE)
individuals <- sapply(fitted_models, function(x) { idx <- regexpr(regexpr_str, x); return(substr(x, idx-3, idx-1)) } )
names(individuals) <- NULL

# get dimension D
fn <- paste0("subsetted_indiv_data/",level,"/",individuals[1],regexpr_str)
load(fn)
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
  load(fn)
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
  distance_mat <- matrix(NA, n_samples_subset*n_indiv, n_samples_subset*n_indiv)
  for(i in 1:(n_indiv*n_samples_subset)) {
    for(j in i:(n_indiv*n_samples_subset)) {
      i_idx <- (i-1)*(P)
      A <- all_samples[,(i_idx+1):(i_idx+(P))]
      j_idx <- (j-1)*(P)
      B <- all_samples[,(j_idx+1):(j_idx+(P))]
      distance_mat[i,j] <- mat_dist(A, B, use_Riemann=use_Riemann)
      distance_mat[j,i] <- distance_mat[i,j]
    }
  }
} else {
  distance_mat <- dist(all_samples)
}

fit <- cmdscale(distance_mat, eig=TRUE, k=3) # k is the number of dim
cat("Eigval 1:",fit$eig[1],"\n")
cat("Eigval 2:",fit$eig[2],"\n")
cat("Eigval 3:",fit$eig[3],"\n")
cat("Eigval 4:",fit$eig[4],"\n")

df <- data.frame(x=fit$points[,1], y=fit$points[,2], z=fit$points[,3], labels=indiv_labels)
save(df, file=paste0("plots/basset/",level,"/",which_measure,"_ordination.RData"))

# centroids are useful for labeling plots
df_centroids <- df %>% group_by(labels) %>% summarise(mean_x=mean(x), mean_y=mean(y), mean_z=mean(z))
save(df_centroids, file=paste0("plots/basset/",level,"/",which_measure,"_ordination_centroids.RData"))

