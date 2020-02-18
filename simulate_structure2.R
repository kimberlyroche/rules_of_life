library(mvtnorm)
library(matrixsampling)

D <- 10
H <- 20
Sigma_arr <- list()
all_samples <- NULL
grp_no <- 3
groups <- sample(3, H, replace=T)
group_corr <- matrix(rnorm(D*grp_no)*0.5, D, grp_no)
for(h in 1:H) {
  etas <- matrix(0, 10, N)
  etas[1,] <- rnorm(N)
  for(i in 2:10) {
    etas[i,] <- etas[1,]*group_corr[i,groups[h]] + rnorm(N)*1
  }
  Sigma_arr[[h]] <- cov(t(etas))
  if(is.null(all_samples)) {
    all_samples <- Sigma_arr[[h]]
  } else {
    all_samples <- cbind(all_samples, Sigma_arr[[h]])
  }
}

d <- Riemann_dist_samples_serial(all_samples, H, 1)
embedding <- cmdscale(d, k=2)
df <- data.frame(x=embedding[,1], y=embedding[,2], label=as.factor(groups))
ggplot(df) +
  geom_point(aes(x=x, y=y, color=label), size=2)










