library(matrixsampling)
library(mvtnorm)
library(driver)

set.seed(100)

T <- 300
D <- 5
Q <- 2 # covariates are seasonal oscillating components

varscaledown <- 0.1

upsilon <- D + 2
Xi <- diag(D)*varscaledown + matrix(varscaledown/10, D, D) - diag(D)*varscaledown/10
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
C.0 <- diag(Q)*varscaledown
M.0 <- matrix(0, Q, D)
Theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1] # mean, Cov(rows), Cov(columns) even though specified is confusing in docs

omega.t <- 2*pi/30 # 10 x wet, 10 x dry
G.t <- matrix(c(cos(omega.t), -sin(omega.t), sin(omega.t), cos(omega.t)), 2, 2)
F.t.T <- matrix(c(1, 0), 1, 2)
W.t <- diag(Q)*varscaledown*0.5
gamma.t <- 1

etas <- matrix(0, D, T)
for(t in 1:T) {
  Theta.t <- G.t%*%Theta.t + rmatrixnormal(1, matrix(0, Q, D), W.t, Sigma)[,,1]
  eta.t.T <- F.t.T%*%Theta.t + rmvnorm(1, rep(0, D), gamma.t*Sigma)
  etas[,t] <- t(eta.t.T)
}

eta.df <- gather_array(etas)
ggplot(eta.df, aes(x=dim_2, y=var, color=as.factor(dim_1))) +
  geom_line()

ys <- matrix(0, D+1, T)
for(t in 1:T) {
  alr.inv <- exp(c(etas[,t], 0))
  alr.inv <- alr.inv/sum(alr.inv)
  ys[,t] <- rmultinom(1, 1000, alr.inv)
}
y.df <- gather_array(ys)
ggplot(y.df, aes(x=dim_2, y=var, fill=as.factor(dim_1))) + 
  geom_bar(position="fill", stat="identity", width=1)

# Kalman filter - learn Sigma and eta.1.T ... eta.T.T
# start with precisely the correct parameterizations
Xi.pred <- Xi
upsilon.pred <- upsilon
G.t.pred <- G.t
F.t.T.pred <- F.t.T
gamma.t.pred <- gamma.t
M.t <- M.0
C.t <- C.0
W.t.pred <- W.t
etas.pred <- matrix(0, D, T)
# we need to keep all these for the smoother!
Ms <- array(0, dim=c(Q, D, T))
Cs <- array(0, dim=c(Q, Q, T))
Rs <- array(0, dim=c(Q, Q, T))
As <- array(0, dim=c(Q, D, T))
Thetas <- array(0, dim=c(Q, D, T))
for(t in 1:T) {
  # prior at t
  A.t <- G.t.pred%*%M.t
  As[,,t] <- A.t
  R.t <- G.t.pred%*%C.t%*%t(G.t.pred) + W.t.pred
  Rs[,,t] <- R.t
  Sigma.pred <- rinvwishart(1, upsilon.pred, Xi.pred)[,,1]
  Theta.t <- rmatrixnormal(1, A.t, R.t, Sigma.pred)[,,1]
  Thetas[,,t] <- Theta.t
  # one-step ahead forecast
  f.t.T = F.t.T.pred%*%A.t
  q.t <- as.numeric(gamma.t.pred + F.t.T.pred%*%R.t%*%t(F.t.T.pred))
  # posterior at t
  e.t.T <- etas[,t] - f.t.T
  etas.pred[,t] <- t(f.t.T)
  S.t <- (R.t%*%t(F.t.T.pred))/q.t
  M.t <- A.t + S.t%*%e.t.T
  Ms[,,t] <- M.t
  C.t <- R.t - q.t*S.t%*%t(S.t)
  Cs[,,t] <- C.t
  upsilon.pred.m.1 <- upsilon.pred
  upsilon.pred <- upsilon.pred + 1
  #Xi.pred <- (1/upsilon.pred)*(upsilon.pred.m.1*Xi.pred + (t(e.t.T)%*%e.t.T)/q.t)
  Xi.pred <- Xi.pred + (t(e.t.T)%*%e.t.T)/q.t
  # this update works better because we need to allow Xi to incrase in scale as upsilon does
  # is this down to a weird parameterization of the inverse Wishart in this R package?
}
Sigma.pred <- rinvwishart(1, upsilon.pred, Xi.pred)[,,1]
Theta.t <- rmatrixnormal(1, M.t, C.t, Sigma.pred)[,,1]
Thetas[,,t] <- Theta.t

cat("Trace of true Xi:         ",sum(diag(Xi)),"\n")
cat("Trace of inferred Xi.pred:",sum(diag(Xi.pred)),"\n")

cat("Trace of true Sigma:    ",sum(diag(Sigma)),"\n")
cat("Trace of inferred Sigma:",sum(diag(Sigma.pred)),"\n")
# the structure looks good but the scale of the inferred Sigma is way too small (???)
# Sigma.pred * upsilon.pred is an almost perfect estimate of Sigma

# plot (one) predicted eta
plot(etas[1,], type="l")
lines(etas.pred[1,], col="blue")

# simulation smoother
Thetas.sim <- array(0, dim=c(Q, D, T))
Thetas.sim[,,T] <- Theta.t
for(t in (T-1):1) {
  Z.t <- Cs[,,t]%*%t(G.t.pred)%*%solve(Rs[,,t])
  M.t.star <- Ms[,,t] + Z.t%*%(Thetas[,,t] - As[,,t])
  C.t.star <- Cs[,,t] - Z.t%*%Rs[,,t]%*%t(Z.t)
  C.t.star <- round(C.t.star, digits=10) # numerical issues, ugh
  Thetas.sim[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.pred)[,,1]
}

# plot components 1 & 2
plot(etas[1,], type="l")
lines(Thetas.sim[1,1,], col="blue")
lines(etas[2,], type="l", lty=2)
lines(Thetas.sim[1,2,], col="red", lty=2)

# Sigma is your residual covariance
# estimate seasonal covariance?
# stack thetas 1:T
big.Theta <- matrix(0, 2*T, D)
for(t in 1:T) {
  idx1 <- (t-1)*2 + 1
  idx2 <- idx1 + 1
  big.Theta[idx1:idx2, 1:D] <- Thetas.sim[,,t]
}
est.cov <- t(big.Theta)%*%big.Theta
est.cov <- est.cov/T
# this is definitely not the right object...
















