# this file repeatedly downsamples simulated data and fits the model, helping us visualize
# the effect of sample number on the quality of posterior estimates

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

library(matrixsampling)

# periodic kernel
#   X is Q x N as in other kernels
#   rho is bandwidth
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
}

# simulate counts
N <- 180
D <- 80
Xi <- matrix(0, D-1, D-1)
Xi[1:12,1:12] <- 0.4
upsilon <- D+20
diag(Xi) <- 1
Sigma <- rinvwishart(1, upsilon+50, (upsilon+50-D)*Xi)[,,1]
X <- matrix(seq(1:N), 1, N)
Theta <- function(X) matrix(0, D-1, ncol(X)) # D must exist in workspace

# square exponential kernel parameters; as in fit_GP
dd_se <- 90
dc <- 0.1 # desired minimum correlation
se_sigma <- 1
rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay

# periodic kernel parameters
period <- 365
per_sigma <- 1
rho_per <- 1

se_weight <- 0.9
per_weight <- 0.1

Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
  per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
  (1e-10)*diag(ncol(X)) # pretty arbitrary

# this scale starts to give realistic looking proportions and counts
Sigma <- Sigma
Lambda <- rmatrixnormal(1, Theta(X), Sigma, Gamma(X))[,,1]
Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N))[,,1]
Pi <- driver::alrInv(t(Eta))
#plot(density(Pi))

Y <- matrix(0, D, N)
for(i in 1:N) {
  Y[,i] <- rmultinom(1, rpois(1, 10000), Pi[i,])
}

par(mfrow=c(1,4))
image(cov2cor(Sigma))

percent_keep <- 0.1
cat("Samples retained:",(N*percent_keep),"\n")
downsampled_idx <- sort(sample(N)[1:round(percent_keep*N)])
downsampled_X <- X[1,downsampled_idx,drop=FALSE]
downsampled_Y <- Y[,downsampled_idx]
temp <- basset(Y=downsampled_Y, X=downsampled_X, Theta=Theta, upsilon=D+2, Xi=(upsilon-D)*diag(D-1), Gamma=Gamma, n_samples=0, ret_mean=TRUE)
png(file.path(relative_path,plot_dir,paste0("downsampled_simulated_correlation_",percent_keep,".png")))
image(cov2cor(temp$Sigma[,,1]))
dev.off()

percent_keep <- 0.2
cat("Samples retained:",(N*percent_keep),"\n")
downsampled_idx <- sort(sample(N)[1:round(percent_keep*N)])
downsampled_X <- X[1,downsampled_idx,drop=FALSE]
downsampled_Y <- Y[,downsampled_idx]
temp <- basset(Y=downsampled_Y, X=downsampled_X, Theta=Theta, upsilon=D+2, Xi=(upsilon-D)*diag(D-1), Gamma=Gamma, n_samples=0, ret_mean=TRUE)
png(file.path(relative_path,plot_dir,paste0("downsampled_simulated_correlation_",percent_keep,".png")))
image(cov2cor(temp$Sigma[,,1]))
dev.off()

percent_keep <- 0.3
cat("Samples retained:",(N*percent_keep),"\n")
downsampled_idx <- sort(sample(N)[1:round(percent_keep*N)])
downsampled_X <- X[1,downsampled_idx,drop=FALSE]
downsampled_Y <- Y[,downsampled_idx]
temp <- basset(Y=downsampled_Y, X=downsampled_X, Theta=Theta, upsilon=D+2, Xi=(upsilon-D)*diag(D-1), Gamma=Gamma, n_samples=0, ret_mean=TRUE)
png(file.path(relative_path,plot_dir,paste0("downsampled_simulated_correlation_",percent_keep,".png")))
image(cov2cor(temp$Sigma[,,1]))
dev.off()





