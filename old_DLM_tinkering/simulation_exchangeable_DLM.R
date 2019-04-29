library(matrixsampling)
# library(mvtnorm)
library(driver)
library(ggplot2)

# =======================================================================================
# SIMULATION - uncollapsed model
# =======================================================================================

set.seed(100)

sim <- 1

T <- 100; D <- 10; Q <- 2

varscale <- 0.01
if(sim > 0) {
  varscale <- log(1.001)
}

delta <- 1 # hack; to get some oscillations, we need large variance somewhere
           # one place to put it is in the initial state

upsilon <- D + 10
# add structure to Sigma through Xi
Xi <- matrix(-0.8, D, D)
Xi[1:(D/2),1:(D/2)] <- 0.8
Xi[(D/2 + 1):D,(D/2 + 1):D] <- 0.8
diag(Xi) <- 1
Xi <- Xi*varscale
if(sim > 0) {
  Xi <- diag(D)*varscale
}
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
C.0 <- diag(Q)*varscale
M.0 <- matrix(0, Q, D)

Theta.t <- rmatrixnormal(1, M.0, C.0, delta*Sigma)[,,1]
if(sim == 1) {
  Theta.t <- cbind(matrix(rnorm(1*(D/2), 0.15, 0.05), 1, (D/2)),
                   matrix(rnorm(1*(D/2), -0.15, 0.05), 1, (D/2)))
  Theta.t <- rbind(Theta.t, -Theta.t)
}
if(sim == 2) {
  Theta.t <- matrix(rnorm(D, 0, 0.025), 1, D)
  Theta.t <- rbind(Theta.t, -Theta.t)
}
Theta.0 <- Theta.t

# want to think about how to make the oscillations themselves noisier

omega.t <- 2*pi/30
G.t <- matrix(c(cos(omega.t), -sin(omega.t), sin(omega.t), cos(omega.t)), 2, 2)
F.t.T <- matrix(c(1, 0), 1, 2)
W.t <- diag(Q)*varscale*0.1
if(sim == 2) {
  W.t <- diag(Q)*varscale*5
}
gamma.t <- 1

etas <- matrix(0, D, T)
Thetas.t <- array(0, dim=c(Q, D, T))
for(t in 1:T) {
  Omega.t <- rmatrixnormal(1, matrix(0, Q, D), W.t, Sigma)[,,1]
  #Theta.t <- G.t%*%Theta.t
  Theta.t <- G.t%*%Theta.t + Omega.t
  Thetas.t[,,t] <- Theta.t
  eta.t.T <- F.t.T%*%Theta.t + rmvnorm(1, rep(0, D), gamma.t*Sigma)
  etas[,t] <- t(eta.t.T)
}

big.Theta <- matrix(0, 2*T, D)
for(t in 1:T) {
  idx1 <- (t-1)*2 + 1
  idx2 <- idx1 + 1
  big.Theta[idx1:idx2, 1:D] <- Thetas.t[,,t]
}
true.theta.cov <- t(big.Theta)%*%big.Theta
true.theta.cov <- true.theta.cov/T
#image(t(true.theta.cov))
png("sim_01_theta_cov_true.png")
heatmap(true.theta.cov, Colv=NA, Rowv=NA)
dev.off()

eta.df <- gather_array(etas[,1:T])
p <- ggplot(eta.df, aes(x=dim_2, y=var, color=as.factor(dim_1))) +
  geom_line()
p
#ggsave("sim_02_eta_oscillating.png", width=8, height=4, units="in", scale=1.5)

ys <- matrix(0, D+1, T)
for(t in 1:T) {
  alr.inv <- exp(c(etas[,t], 0))
  alr.inv <- alr.inv/sum(alr.inv)
  ys[,t] <- rmultinom(1, 1000, alr.inv)
}
y.df <- gather_array(ys[,1:T])
p <- ggplot(y.df, aes(x=dim_2, y=var, fill=as.factor(dim_1))) + 
  geom_bar(position="fill", stat="identity", width=1)
p
#ggsave("sim_03_y_oscillating.png", width=8, height=4, units="in", scale=1.5)

sum(diag(W.t))
sum(diag(Sigma))

# =======================================================================================
# SIMULATION - collapsed model
# =======================================================================================

B <- matrix(0, D, T)
for(t in 1:T) {
  alpha.t <- F.t.T
  for(j in t:1) {
    alpha.t <- alpha.t%*%G.t
  }
  #alpha.t <- t(alpha.t%*%M.0)
  alpha.t <- t(alpha.t%*%Theta.0)
  B[,t] <- alpha.t
}

K <- Xi

A <- matrix(0, T, T)
# because G.t is fixed, I'm just power iterating these
for(t in 1:T) {
  for(tk in t:T) {
    k <- tk - t
    if(t == tk) {
      sum_val <- W.t
      for(ell in t:2) {
        sum_val <- sum_val + G.t**abs(t - ell)%*%W.t%*%(t(G.t)**abs(ell - t))
      }
      sum_val <- sum_val + delta*G.t**abs(t - 1)%*%C.0%*%(t(G.t)**abs(t - 1))
      sum_val <- F.t.T%*%sum_val%*%t(F.t.T)
      A[t, tk] <- gamma.t + sum_val
    } else {
      sum_val <- G.t**abs(t - (t - k + 1))%*%W.t
      for(ell in (t-k):2) {
        sum_val <- sum_val + G.t**abs(t - ell)%*%W.t%*%(t(G.t)**abs(ell - (t-k)))
      }
      sum_val <- sum_val + delta*G.t**abs(t - 1)%*%C.0%*%(t(G.t)**abs(1 - (t-k)))
      A[t, tk] <- F.t.T%*%sum_val%*%t(F.t.T)
    }
  }
}
for(t in 1:T) {
  for(tk in 1:(t-1)) {
    A[t, tk] <- A[tk, t]
  }
}

# this sucks! A is not always positive semi-definite; depends on gamma.t, G.t, W.t relative magnitudes
eigen(A)$values

etas.collapsed <- rmatrixt(1, upsilon, B, K, A)[,,1]

eta.collapsed.df <- gather_array(etas.collapsed[,1:T])
p <- ggplot(eta.collapsed.df, aes(x=dim_2, y=var, color=as.factor(dim_1))) +
  geom_line()
p
#ggsave("sim_02b_eta_oscillating.png", width=8, height=4, units="in", scale=1.5)

ys.collapsed <- matrix(0, D+1, T)
for(t in 1:T) {
  alr.inv <- exp(c(etas.collapsed[,t], 0))
  alr.inv <- alr.inv/sum(alr.inv)
  ys.collapsed[,t] <- rmultinom(1, 1000, alr.inv)
}
y.collapsed.df <- gather_array(ys.collapsed[,1:T])
p <- ggplot(y.collapsed.df, aes(x=dim_2, y=var, fill=as.factor(dim_1))) + 
  geom_bar(position="fill", stat="identity", width=1)
p
#ggsave("sim_03b_y_oscillating.png", width=8, height=4, units="in", scale=1.5)

# =======================================================================================
# PREDICTION
# =======================================================================================

# Kalman filter - learn Sigma and eta.1.T ... eta.T.T
# start with precisely the correct parameterizations
Xi.pred <- Xi
upsilon.pred <- upsilon
G.t.pred <- G.t
F.t.T.pred <- F.t.T
gamma.t.pred <- gamma.t
M.t <- M.0
C.t <- diag(Q)*10
W.t.pred <- array(rep(W.t, T), dim=c(Q, Q, T))
for(i in 1:30) {
  W.t.pred[,,i] <- diag(Q)
}

etas.pred <- matrix(0, D, T)
# we need to keep all these for the smoother!
Ms <- array(0, dim=c(Q, D, T))
Cs <- array(0, dim=c(Q, Q, T))
Rs <- array(0, dim=c(Q, Q, T))
As <- array(0, dim=c(Q, D, T))
Thetas <- array(0, dim=c(Q, D, T))

# the trouble is this is very sensitive to *where* in the oscillation a given component starts
# if it's initialized in the wrong spot (i.e. out of phase), all the error gets pushed into Sigma
# (the residual)

# making the state transition variance large in initial states helps enormously

for(t in 1:T) {
  # prior at t
  A.t <- G.t.pred%*%M.t
  As[,,t] <- A.t
  R.t <- round(G.t.pred%*%C.t%*%t(G.t.pred) + W.t.pred[,,t], digits=10) # numerical issues, need a better hack
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

png("sim_04_Sigma_true.png")
heatmap(Sigma, Colv=NA, Rowv=NA)
dev.off()

png("sim_05_Sigma_pred.png")
heatmap(Sigma.pred, Colv=NA, Rowv=NA)
dev.off()

Thetas.sim <- array(0, dim=c(Q, D, T))
Thetas.sim[,,T] <- Theta.t
for(t in (T-1):1) {
  Z.t <- Cs[,,t]%*%t(G.t.pred)%*%solve(Rs[,,t])
  M.t.star <- Ms[,,t] + Z.t%*%(Thetas[,,t] - As[,,t])
  C.t.star <- Cs[,,t] - Z.t%*%Rs[,,t]%*%t(Z.t)
  C.t.star <- round(C.t.star, digits=10) # numerical issues, ugh
  Thetas.sim[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.pred)[,,1]
}

d <- 1
temp.df <- data.frame(rbind(data.frame(t=1:T, val=etas[d,1:T], thing=rep("eta", T)),
                            data.frame(t=1:T, val=Thetas.sim[1,d,1:T], thing=rep("Theta.sim", T))))
p <- ggplot(temp.df, aes(x=t, y=val, color=thing)) +
  geom_line()
ggsave(paste("sim_06_eta-theta_pair.png", sep=""), width=8, height=4, units="in", scale=1.5)

big.Theta <- matrix(0, 2*T, D)
for(t in 1:T) {
  idx1 <- (t-1)*2 + 1
  idx2 <- idx1 + 1
  big.Theta[idx1:idx2, 1:D] <- Thetas.sim[,,t]
}
est.cov <- t(big.Theta)%*%big.Theta
est.cov <- est.cov/T
png("sim_07_theta_cov_est.png")
heatmap(est.cov, Colv=NA, Rowv=NA)
dev.off()

temp.df <- rbind(data.frame(x=1:T, var=etas[1,], thing="eta1"), data.frame(x=1:T, var=etas[2,], thing="eta2"))
p <- ggplot(temp.df, aes(x=x, y=var, color=thing)) +
  geom_line()
ggsave(paste("sim_08_high_pcorr_etas.png", sep=""), width=8, height=4, units="in", scale=1.5)

temp.df <- rbind(data.frame(x=1:T, var=etas[1,], thing="eta1"), data.frame(x=1:T, var=etas[6,], thing="eta2"))
p <- ggplot(temp.df, aes(x=x, y=var, color=thing)) +
  geom_line()
ggsave(paste("sim_09_high_ncorr_etas.png", sep=""), width=8, height=4, units="in", scale=1.5)

# windowed correlation
#idx1 <- 1
#idx2 <- 2
#window <- 30
#cor_vec <- numeric(T-window)
#for(t in 1:(T-window)) {
#  rho <- stats::cor(c(Thetas.sim[,idx1,t:(t+window)]), c(Thetas.sim[,idx2,t:(t+window)]))
#  cor_vec[t] <- rho
#}
#plot(cor_vec, type="l", ylim=c(-1,1))
























