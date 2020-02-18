library(matrixsampling)
library(stray)
library(driver)
library(Rcpp)
library(ggplot2)

sourceCpp("include/cpp/Riemann_dist.cpp")

# ===============================================================================================
#   SIMULATION 1: vanilla simulation directly from model
# ===============================================================================================

D <- 50
H <- 5
n_samples <- 10
Y_arr <- list()
X_arr <- list()
Sigma_arr <- list()
fit_arr <- list()

# priors
upsilon <- D + 10
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

for(h in 1:H) {
  cat("Simulating individual",h,"\n")
  Sigma_arr[[h]] <- rinvwishart(1, upsilon, Xi)[,,1]
  N <- rpois(1, 75)
  X <- matrix(sort(sample(1:1000)[1:N]), 1, N)
  Theta <- function(X) matrix(0, D-1, ncol(X)) # not particularly realistic
  Gamma <- function(X) SE(X, sigma=1, rho=10)
  Lambda <- rmatrixnormal(1, Theta(X), Sigma_arr[[h]], Gamma(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma_arr[[h]], diag(N))[,,1]
  pi <- t(alrInv(t(Eta))) # D x N
  X_arr[[h]] <- X
  Y_arr[[h]] <- matrix(0, D, N)
  for(i in 1:N) {
    Y_arr[[h]][,i] <- rmultinom(1, rpois(1, 10000), pi[,i])
  }
  fit_arr[[h]] <- basset(Y_arr[[h]], X_arr[[h]], upsilon, Theta, Gamma, Xi, n_samples=n_samples)
}

par(mfrow=c(H,2))
for(h in 1:H) {
  image(Sigma_arr[[h]])
  image(apply(fit_arr[[h]]$Sigma, c(1,2), mean))
}

all_samples <- matrix(NA, D-1, H*(D-1)*n_samples)
labels <- c()
for(h in 1:H) {
  host_offset <- (h-1)*(D-1)*n_samples
  labels <- c(labels, rep(h, n_samples))
  for(j in 1:n_samples) {
    sample_offset <- (j-1)*(D-1)
    all_samples[,(host_offset+sample_offset+1):(host_offset+sample_offset+D-1)] <- fit_arr[[h]]$Sigma[,,j]
  }
}

d <- Riemann_dist_samples_serial(all_samples, H, n_samples)
embedding <- cmdscale(d, k=2)
df <- data.frame(x=embedding[,1], y=embedding[,2], label=as.factor(labels))
ggplot(df) +
  geom_point(aes(x=x, y=y, color=label))

# ===============================================================================================
#   SIMULATION 2: different Sigma.1 (moving baseline dynamics) and Sigma.2 (residual)
#                 (1) some correlation at moving baseline level (lambda level)
#                 (2) some correlation in deviations from baseline (eta level)
#                 (3) same Sigma (w/ correlation)
# ===============================================================================================

set.seed(1)

D <- 50
H <- 4
n_samples <- 10
Y_arr <- list()
X_arr <- list()
Sigma_arr.1 <- list()
Sigma_arr.2 <- list()
fit_arr <- list()

which_sim <- 2

# priors
upsilon <- D + 10

GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi_corr <- GG%*%(diag(D)*0.4 + 0.6)%*%t(GG) # take diag as covariance over log abundances
Xi_corr <- Xi_corr*(upsilon-D-1)

for(h in 1:H) {
  cat("Simulating individual",h,"\n")
  Sigma_arr.1[[h]] <- rinvwishart(1, upsilon, Xi)[,,1]
  Sigma_arr.2[[h]] <- rinvwishart(1, upsilon, Xi_corr)[,,1]
  N <- rpois(1, 75)
  X <- matrix(sort(sample(1:1000)[1:N]), 1, N)
  Theta <- function(X) matrix(0, D-1, ncol(X)) # not particularly realistic
  Gamma <- function(X) SE(X, sigma=1, rho=10)
  if(which_sim == 1) {
    Lambda <- rmatrixnormal(1, Theta(X), Sigma_arr.2[[h]], Gamma(X))[,,1]
    #Eta <- rmatrixnormal(1, Lambda, Sigma_arr.1[[h]], diag(N))[,,1]
    Eta <- rmatrixnormal(1, Lambda, diag(D-1), diag(N))[,,1]
  }
  if(which_sim == 2) {
    #Lambda <- rmatrixnormal(1, Theta(X), Sigma_arr.1[[h]], Gamma(X))[,,1]
    Lambda <- rmatrixnormal(1, Theta(X), diag(D-1), Gamma(X))[,,1]
    Eta <- rmatrixnormal(1, Lambda, Sigma_arr.2[[h]], diag(N))[,,1]
  }
  if(which_sim == 3) {
    Lambda <- rmatrixnormal(1, Theta(X), Sigma_arr.2[[h]], Gamma(X))[,,1]
    Eta <- rmatrixnormal(1, Lambda, Sigma_arr.2[[h]], diag(N))[,,1]
  }
  pi <- t(alrInv(t(Eta))) # D x N
  X_arr[[h]] <- X
  Y_arr[[h]] <- matrix(0, D, N)
  for(i in 1:N) {
    Y_arr[[h]][,i] <- rmultinom(1, rpois(1, 10000), pi[,i])
  }
  fit_arr[[h]] <- basset(Y_arr[[h]], X_arr[[h]], upsilon, Theta, Gamma, Xi, n_samples=n_samples)
}

par(mfrow=c(H,2))
for(h in 1:H) {
  image(Sigma_arr.2[[h]])
  image(apply(fit_arr[[h]]$Sigma, c(1,2), mean))
}

all_samples <- matrix(NA, D-1, H*(D-1)*n_samples)
labels <- c()
for(h in 1:H) {
  host_offset <- (h-1)*(D-1)*n_samples
  labels <- c(labels, rep(h, n_samples))
  for(j in 1:n_samples) {
    sample_offset <- (j-1)*(D-1)
    all_samples[,(host_offset+sample_offset+1):(host_offset+sample_offset+D-1)] <- fit_arr[[h]]$Sigma[,,j]
  }
}

d <- Riemann_dist_samples_serial(all_samples, H, n_samples)
embedding <- cmdscale(d, k=2)
df <- data.frame(x=embedding[,1], y=embedding[,2], label=as.factor(labels))
ggplot(df) +
  geom_point(aes(x=x, y=y, color=label))





