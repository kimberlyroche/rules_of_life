library(matrixsampling)
library(mvtnorm)
library(driver)

set.seed(102)

T <- 300
D <- 10
Q <- 2 # covariates are seasonal oscillating components

varscaledown <- 0.1

upsilon <- D + 2
Xi <- diag(D)*varscaledown + matrix(varscaledown/10, D, D) - diag(D)*varscaledown/10

# structure in Sigma - does this come through? is it interpretable?
# basically: we want to impose structure between log ratios AND recover it

# the fact that you can get ~0.5 correlation out of random interactions should be a point for caution though...

#Xi <- matrix(-varscaledown/10, D, D)
#Xi[1:5,1:5] <- varscaledown/10
#Xi[6:10,6:10] <- varscaledown/10
#diag(Xi) <- varscaledown

Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
C.0 <- diag(Q)*varscaledown
M.0 <- matrix(0, Q, D)
Theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1] # mean, Cov(rows), Cov(columns) even though specified is confusing in docs

omega.t <- 2*pi/30
G.t <- matrix(c(cos(omega.t), -sin(omega.t), sin(omega.t), cos(omega.t)), 2, 2)
F.t.T <- matrix(c(1, 0), 1, 2)
W.t <- diag(Q)*varscaledown*0.5 # reducing the scale of W.t means things stay almost perfectly in phase
                                # increasing the scale of W.t allows things to get out of phase
                                # BUT the length of period is pretty rigid here, do we want to change that?
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
y.df <- gather_array(ys[,100:T])
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
  Xi.pred <- Xi.pred + (t(e.t.T)%*%e.t.T)/q.t
}
Sigma.pred <- rinvwishart(1, upsilon.pred, Xi.pred)[,,1]
Sigma.pred.cor <- cov2cor(Sigma.pred)
Theta.t <- rmatrixnormal(1, M.t, C.t, Sigma.pred)[,,1]
Thetas[,,t] <- Theta.t

cat("Trace of true Sigma:    ",sum(diag(Sigma)),"\n")
cat("Trace of inferred Sigma:",sum(diag(Sigma.pred)),"\n")

# plot predicted etas
for(d in 1:D) {
  plot(etas[d,100:T], type="l")
  lines(etas.pred[d,100:T], col="blue")
}

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

for(d in 1:D) {
  plot(etas[d,], type="l")
  lines(Thetas.sim[1,d,], col="blue")
}

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
est.cor <- cov2cor(est.cov)

print(round(est.cor, digits=5))
image(est.cor)
print(round(Sigma.pred.cor, digits=5))
image(Sigma.pred.cor)

# estimated seasonal correlation pairs
# 1 vs. 2 have approx. +0.5 correlation
plot(etas[1,], type="l")
lines(etas[4,], col="red")

# 1 vs. 4 have approx. 0 correlation
plot(etas[1,], type="l")
lines(etas[6,], col="red")

# 1 vs. 5 have approx -0.5 correlation
plot(etas[2,], type="l")
lines(etas[5,], col="red")

# residual correlation pairs are super hard to see/interpret

# if we assume these (seasonal) covariances are what we want: how variable are these over the time course?
# 6 and 9 have very high correlation
# should generally be the case that the only way to achieve that will be by always being in phase!

# what does windowed correlation look like?
idx1 <- 6
idx2 <- 9
window <- 30
cor_vec <- numeric(T-window)
for(t in 1:(T-window)) {
  rho <- stats::cor(c(Thetas.sim[,idx1,t:(t+window)]), c(Thetas.sim[,idx2,t:(t+window)]))
  cor_vec[t] <- rho
}
plot(cor_vec, type="l", ylim=c(-1,1))
# these can vary a lot, as they should, because we're not building in any particular covariance (?)















