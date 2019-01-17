library(mvtnorm)
library(MCMCpack)
library(LaplacesDemon)

# marginal matrix-T (two component)
logd_mat_t_two <- function(s, vc, data, N, P, upsilon) {
  K <- diag(P)
  A <- round((exp(s[1])*vc[[1]] + exp(s[2])*vc[[2]]), digits=10)
  # garbage workaround for numeric precision problems
  samples <- dim(data)[3]
  d <- 0
  for(i in 1:samples) {
    d <- d - (P/2)*log(det(A)) - ((upsilon+N+P-1)/2)*log(det(diag(P) + (1/upsilon)*solve(K)%*%(data[,,i])%*%solve(A)%*%t(data[,,i])))
  }
  return(-d)
}

N <- 300
P <- 20
samples <- 1000
upsilon <- P + 2
M <- matrix(0, P, N)

# variance component 1
sigma_true.1 <- 0.2
component.1 <- diag(N)
for(i in 1:(N-1)) {
  component.1[i,i+1] <- 0.33
  component.1[i+1,i] <- 0.33
}

# variance component 2
sigma_true.2 <- 0.8
component.2 <- matrix(-0.2, nrow=N, ncol=N)
component.2[1:150,1:150] <- 0.2
component.2[151:300,151:300] <- 0.2
for(i in 1:N) {
  component.2[i,i] <- 1
}

sample_data <- array(0, dim=c(P, N, samples))
V <- round((sigma_true.1*component.1 + sigma_true.2*component.2), digits=10)
for(i in 1:samples) {
  U <- round(riwish(upsilon, diag(P)), digits=5)
  sample_data[,,i] <- rmatrixnorm(M, U, V)
}

res <- optim(par=c(runif(1), runif(1)),
             fn=logd_mat_t_two,
             vc=list(component.1, component.2),
             data=sample_data, N=N, P=P, upsilon=upsilon)
cat("Actual:",sigma_true.1,", Estimated:",exp(res$par[1]),"\n")
cat("Actual:",sigma_true.2,", Estimated:",exp(res$par[2]),"\n")
total_wgt_actual <- sigma_true.1 + sigma_true.2
total_wgt <- exp(res$par[1]) + exp(res$par[2])
cat("Actual proportion (1):",(sigma_true.1/total_wgt_actual),
    ", Estimated proportion (1):",(exp(res$par[1]))/total_wgt,"\n")
cat("Actual proportion (2):",(sigma_true.2/total_wgt_actual),
    ", Estimated proportion (2):",(exp(res$par[2]))/total_wgt,"\n")










