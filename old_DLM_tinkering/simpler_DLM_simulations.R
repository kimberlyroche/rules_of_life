library(mvtnorm)
library(matrixsampling)

# univariate version
# single harmonic

F <- matrix(c(1, 0), 1, 2)
theta.t <- matrix(c(1, 1), 2, 1)
V <- 0.25 # sd
omega <- 2*pi/10
G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
W <- diag(2)*0.25
T <- 100

ys <- numeric(T)
for(t in 1:T) {
  theta.t <- G%*%theta.t + rnorm(1, 0, W)
  ys[t] <- F%*%theta.t + rnorm(1, 0, V)
}
plot(ys, type="l")

# multivariate (p=3) version
# single harmonic

theta.t <- matrix(rnorm(2), 2, 3)
V <- diag(3)*0.1 # sd
omega <- 2*pi/10
G <- diag(2)
W <- diag(2)*0.05
T <- 100

ys <- matrix(0, T, 3)
Fts <- matrix(0, T, 2)
for(t in 1:T) {
  theta.t <- G%*%theta.t + rmatrixnormal(1, matrix(0, 2, 3), W, V)[,,1]
  F.t <- matrix(c(cos(omega*t), sin(omega*t)), 1, 2)
  Fts[t,] <- F.t
  Fts[t] <- F.t%*%matrix(1, 2, 1)
  ys[t,] <- F.t%*%theta.t + rmvnorm(1, rep(0, 3), V)
}
plot(ys[,1], type="l")
lines(ys[,2], col="blue")
lines(ys[,3], col="red")

# try covariance calculation on simple version

C.0 <- W
covmat <- matrix(0, T, T)
gamma.t <- 1
for(i in 1:T) {
  for(j in 1:T) {
    if(j >= i) {
      if(i == j) {
        temp <- W
        for(ell in j:2) {
          power <- length(j:ell)
          temp <- temp + (G**power)%*%W%*%(t(G)**power)
        }
        power <- length(j:1)
        temp <- temp + (G**power)%*%C.0%*%(t(G)**power)
        covmat[i,j] <- gamma.t + Fts[j,]%*%temp%*%matrix(Fts[j,], 2, 1)
      } else {
        power <- length(j:(i+1))
        temp <- (G**power)%*%W
        for(ell in i:2) {
          power <- length(j:ell)
          power2 <- length(ell:i)
          temp <- temp + (G**power)%*%W%*%(t(G)**power2)
        }
        power <- length(j:1)
        power2 <- length(1:i)
        temp <- temp + (G**power)%*%C.0%*%(t(G)**power2)
        covmat[i,j] <- Fts[j,]%*%temp%*%matrix(Fts[j,], 2, 1)
      }
    } else {
      covmat[j,i] <- covmat[i,j]
    }
  }
}









