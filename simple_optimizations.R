library(mvtnorm)
library(MCMCpack)
library(LaplacesDemon)

# univariate normal

sigma_true <- 2
N <- 100
data <- rnorm(100, 0, sigma_true)

logd <- function(scale) {
  sigma <- scale
  # assume mean zero
  d <- -(N/2)*log(2*pi*sigma^2) - sum(unlist(lapply(data, function(x, y) { (1/2)*(x/y)^2 }, y=sigma )))
  return(-d)
}

avg_est <- 0
for(i in 1:10) {
  res <- optim(par=c(runif(1)), fn=logd, method="L-BFGS-B")
  avg_est <- avg_est + res$par
}
cat("Actual:",sigma_true,", avg est:",avg_est/10,"\n")

# multivariate normal

D <- 5
Sigma_true <- 2*diag(D)
data <- t(rmvnorm(N, mean=rep(0,D), sigma=Sigma_true))

logd_mult <- function(scale) {
  Sigma <- scale*diag(D)
  # assume mean zero
  #d <- -(N/2)*log(det(2*pi*Sigma)) - (1/2)*t(data)%*%solve(Sigma)%*%data
  d <- -(N/2)*log(det(2*pi*Sigma)) - sum(apply(data, 2, function(x, y) { (1/2)*t(x)%*%solve(y)%*%x }, y=Sigma))
  return(-d)
}

avg_est <- 0
for(i in 1:10) {
  res <- optim(par=c(runif(1)), fn=logd_mult, method="L-BFGS-B")
  avg_est <- avg_est + res$par
}
cat("Actual:\n")
print(Sigma_true)
cat("Avg est:",(avg_est/10),"\n")

# matrix normal

P <- 20
N <- 30
it <- 500
zero_mean <- matrix(0, P, N)
sigma_true.1 <- 0.25
component.1 <- diag(N)
sigma_true.2 <- 2
component.2 <- matrix(0, nrow=N, ncol=N)
component.2[1:15,1:15] <- 0.8
component.2[16:30,16:30] <- 0.8
component.2 <- component.2 + diag(N)*0.1
image(component.2)
sample_data <- array(0, dim=c(P, N, it))
sample_covar <- sigma_true.1*component.1 + sigma_true.2*component.2
image(sample_covar)
for(i in 1:it) {
  sample_data[,,i] <- rmatrixnorm(zero_mean, diag(P), sample_covar)
}
image(sample_data[,,1])

logd_mat_norm <- function(s) {
  s1 <- s[1]
  s2 <- s[2]
  A <- round((exp(s1)*component.1 + exp(s2)*component.2), digits=10)
  d <- 0
  for(i in 1:it) {
    d <- d - dmatrixnorm(sample_data[,,i], rep(0, P), diag(P), A, log=TRUE)
  }
  return(d)
}

res <- optim(par=c(runif(1), runif(1)), fn=logd_mat_norm, method="L-BFGS-B")
cat("Actual:",sigma_true.1,",",sigma_true.2," avg est:",exp(res$par[1]),",",exp(res$par[2]),"\n")

# marginal matrix T - one component

it <- 200
P <- 15
N <- 20
upsilon <- P + 100
zero_mean <- matrix(0, P, N)
sigma_true.1 <- 2
component.1 <- diag(N)
sample_data <- array(0, dim=c(P, N, it))
sample_covar.2 <- round((sigma_true.1*component.1), digits=10)
for(i in 1:it) {
  sample_covar.1 <- round(riwish(upsilon, diag(P)), digits=10)
  sample_data[,,i] <- rmatrixnorm(zero_mean, sample_covar.1, sample_covar.2)
}

logd_mat_t <- function(s) {
  s1 <- s[1]
  K <- diag(P)
  A <- round((exp(s1)*component.1), digits=10)
  d <- 0
  for(j in 1:it) {
    d <- d + (P/2)*log(det(A)) +
      ((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(sample_data[,,i])%*%solve(A)%*%t(sample_data[,,i])))
  }
  return(d)
}

K <- diag(P)
res <- optim(par=c(runif(1)), fn=logd_mat_t, method="L-BFGS-B")
cat("Actual:",sigma_true.1,", avg est:",exp(res$par[1]),"\n")

# Q1: does running this for 100,000 iterations get super close to the right scale?

# marginal matrix T - two components

it <- 500
P <- 15
N <- 20
upsilon <- P + 100
zero_mean <- matrix(0, P, N)
sigma_true.1 <- 0.25
component.1 <- diag(N)
sigma_true.2 <- 2
component.2 <- matrix(0, nrow=N, ncol=N)
component.2[1:10,1:10] <- 0.8
component.2[11:20,11:20] <- 0.8
component.2 <- component.2 + diag(N)*0.1
sample_data <- array(0, dim=c(P, N, it))
sample_covar.2 <- round((sigma_true.1*component.1 + sigma_true.2*component.2), digits=10)
for(i in 1:it) {
  sample_covar.1 <- round(riwish(upsilon, diag(P)), digits=10)
  sample_data[,,i] <- rmatrixnorm(zero_mean, sample_covar.1, sample_covar.2)
}

logd_mat_t <- function(s) {
  s1 <- s[1]
  s2 <- s[2]
  K <- diag(P)
  A <- round((exp(s1)*component.1 + exp(s2)*component.2), digits=10)
  d <- 0
  for(j in 1:it) {
    d <- d + (P/2)*log(det(A)) +
      ((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(sample_data[,,i])%*%solve(A)%*%t(sample_data[,,i])))
  }
  return(d)
}

K <- diag(P)
res <- optim(par=c(runif(1), runif(1)), fn=logd_mat_t, method="L-BFGS-B")
cat("Actual:",sigma_true.1,",",sigma_true.2," avg est:",exp(res$par[1]),",",exp(res$par[2]),"\n")













