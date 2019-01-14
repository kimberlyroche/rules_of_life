library(mvtnorm)
library(MCMCpack)
library(LaplacesDemon)

# DENSITIES

# log univariate normal density
logd_uni_normal <- function(s, data, N) {
  # assume mean zero
  d <- 0
  for(i in 1:N) {
    d <- d + dnorm(data[i], 0, s[1], log=TRUE)
  }
  return(-d)
}

# log multivariate normal density
logd_mult_normal <- function(s, data, N, D) {
  Sigma <- s[1]*diag(D)
  d <- -(N/2)*log(det(2*pi*Sigma)) - sum(apply(data, 2, function(x, y) { (1/2)*t(x)%*%solve(y)%*%x }, y=Sigma))
  return(-d)
}

# matrix normal
logd_mat_norm <- function(s, vc, data, P) {
  s1 <- s[1]
  s2 <- s[2]
  A <- round((exp(s1)*vc[[1]] + exp(s2)*vc[[2]]), digits=10)
  d <- 0
  samples <- dim(data)[3]
  for(i in 1:samples) {
    d <- d + dmatrixnorm(data[,,i], rep(0, P), diag(P), A, log=TRUE)
  }
  return(-d)
}

# marginal matrix-T (one component)
logd_mat_t_one <- function(s, vc, data, N, P, upsilon) {
  K <- diag(P)
  A <- round((exp(s[1])*vc[[1]]), digits=10)
  samples <- dim(data)[3]
  d <- 0
  for(i in 1:samples) {
    d <- d - (P/2)*log(det(A)) -
      ((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(data[,,i])%*%solve(A)%*%t(data[,,i])))
  }
  return(-d)
}

# marginal matrix-T (two component)
logd_mat_t_two <- function(s, vc, data, N, P, upsilon) {
  K <- diag(P)
  A <- round((exp(s[1])*vc[[1]] + exp(s[2])*vc[[2]]), digits=5)
  samples <- dim(data)[3]
  d <- 0
  for(i in 1:samples) {
    d <- d - (P/2)*log(det(A)) -
      ((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(data[,,i])%*%solve(A)%*%t(data[,,i])))
  }
  return(-d)
}

