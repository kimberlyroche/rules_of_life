# this file simulates adding noise to a PD matrix and track changes in Riemannian distance to better understand
# "small" vs. "large" distances in this space; here we add noise by sampling from an inverse Wishart with
# smaller upsilon (larger degrees of freedom?)
#
# this is replaced by `simulate_Riemannian_ATA.R` because it's very hard to interpret the change in distance
# relative to the non-linear change in average element-wise variance created by a unit change in upsilon
#
# the derivative of d(Var)/d(upsilon) (calculated below) was an attempt to warp the x-axis of a plot of
# upsilon vs. Riemannian distance to make this interpretable but it's messy (and probably incorrect?)

relative_path <- ".."

library(Rcpp)
library(matrixsampling)
library(driver)
library(ggplot2)
library(gridExtra)

sourceCpp(file.path(relative_path,"include/cpp/Riemann_dist.cpp"))

# derivative of variance in an arbitrary element of a matrix ~ IW relative to
# a change in upsilon; we'll use these functions to warp the x-axis of the
# noise x average distance plot such that we have something like a constant
# increase in noise; in reality, noise is a non-linear function of upsilon
num_deriv <- function(Xi, i, j) {
  return(Xi[i,j]^2 + Xi[i,i]*Xi[j,j])
}
denom_deriv <- function(upsilon, D) {
  P <- D - 1
  (upsilon - P - 1)*(4*(upsilon^2) + (-8*P - 11)*upsilon + 4*(P^2) + 11*P + 3)
}
IW_deriv <- function(Xi, upsilon, D, i, j) {
  P <- D - 1
  A <- (upsilon - P + 1)*(Xi[i,j]^2) + (upsilon - P - 1)*Xi[i,i]*Xi[j,j]
  A_deriv <- num_deriv(Xi, i, j)
  B <- (upsilon - P)*((upsilon - P - 1)^2)*(upsilon - P - 3)
  B_deriv <- denom_deriv(upsilon, D)
  return((B*A_deriv - A*B_deriv)/(B^2))
}

D <- 20 # dimension of covariance matrix (D x D)
N <- 50 # number of noisy additions
# sample a baseline covariance matrix
samples <- NULL
upsilon <- D + 10
baseline <- rinvwishart(1, upsilon, diag(D)*(upsilon - D - 1))[,,1] # baseline
# vector to store distances in
mean_distances <- numeric(N-1)
lower_distances <- numeric(N-1)
upper_distances <- numeric(N-1)
mean_distances_Euclidean <- numeric(N-1)
lower_distances_Euclidean <- numeric(N-1)
upper_distances_Euclidean <- numeric(N-1)
steps <- c()
step_sizes <- c()
for(i in 1:(N-1)) {
  # successively add more positive (semi) definite noise to a baseline covariance matrix
  # by moving away from the mean
  upsilon <- D + 2 + (N - i) # note: min upsilon is D + 2
  cat("Upsilon:",upsilon,"\n")
  steps <- c(steps, upsilon)
  step_sizes <- c(step_sizes, -IW_deriv(baseline, upsilon, D, 1, 2)) # negative because we're interested change relative
                                                                     # to a negative unit increase in upsilon
  dvec <- c()
  dvec_RMSE <- c()
  for(j in 1:100) {
    new_sample <- rinvwishart(1, upsilon, baseline*(upsilon - D - 1))[,,1]
    if(is.null(samples)) {
      samples <- new_sample
    } else {
      samples <- cbind(samples, new_sample)
    }
    # calculate distance of this noisy matrix from the original
    dvec <- c(dvec, Riemann_dist_pair(baseline, new_sample))
    dvec_RMSE <- c(dvec_RMSE, sqrt(mean((c(baseline) - c(new_sample))^2)))
    # check PD requirement
    if(min(eigen(new_sample)$values) <= 0) {
      cat("Not positive definite!\n")
    }
  }
  mean_distances[i] <- mean(dvec)
  lower_distances[i] <- sort(dvec)[10]
  upper_distances[i] <- sort(dvec)[90]
  mean_distances_Euclidean[i] <- mean(dvec_RMSE)
  lower_distances_Euclidean[i] <- sort(dvec_RMSE)[10]
  upper_distances_Euclidean[i] <- sort(dvec_RMSE)[90]
}
# normalize steps
step_sizes <- step_sizes - min(step_sizes)
step_sizes <- step_sizes*(1/max(step_sizes))

# replacesd 1:(N-1) with step_sizes
# distances_plot <- data.frame(sample=1:(N-1), mean=mean_distances, lower=lower_distances, upper=upper_distances, type="Riemannian")
# distances_plot <- rbind(distances_plot, data.frame(sample=1:(N-1), mean=mean_distances_Euclidean, lower=lower_distances_Euclidean, upper=upper_distances_Euclidean, type="RMSE"))

# replacesd 1:(N-1) with step_sizes
distances_plot <- data.frame(sample=step_sizes, mean=mean_distances, lower=lower_distances, upper=upper_distances, type="Riemannian")
distances_plot <- rbind(distances_plot, data.frame(sample=step_sizes, mean=mean_distances_Euclidean, lower=lower_distances_Euclidean, upper=upper_distances_Euclidean, type="RMSE"))

# plot distances
ggplot(distances_plot) +
  geom_path(aes(x=sample, y=mean), size=1) +
  geom_ribbon(aes(x=sample, ymin=lower, ymax=upper, alpha=0.33)) +
  facet_wrap(. ~ type, scales = "free_y") +
  xlab("noise addition step")






