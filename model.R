source("include.R")

load("filtered_90.RData")
cat(ntaxa(f2),"x",nsamples(f2),"\n")

DUI_counts <- subset_samples(f2, sname=="DUI")
save(DUI_counts, file="DUI_filtered_counts.RData")

quit()

DUI_log_ratios <- subset_samples(clr_f2, sname=="DUI")

ilr_counts <- apply(otu_table(f2)@.Data+0.65, 1, ilr)
dim(ilr_counts)
rm(f2)

# ilr_counts is now taxa (rows) x samples (columns)
mu <- apply(ilr_counts, 1, mean)


matT_log_density <- function(sigma) {
  upsilon <- P+2
  A <- sigma[1]*diag(N) + sigma[2]*diag(N)*0.5
  # rem redundant parts
  d <- log(gamma((upsilon+N+P-1)/2)) -
       log(gamma((upsilon+P-1)/2)) -
       ((N*P)/2)*log(pi) -
       (N/2)*log(det(K)) -
       (P/2)*log(det(A)) -
       ((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(eta-mu)%*%solve(A)%*%t(eta-mu)))
  return(-d)
}

N <- 30
P <- 27
eta <- ilr_counts[,1:N]
mu <- apply(eta, 1, mean)
K <- cov(t(eta))

#res <- optim(c(0.25, 0.5), matT_log_density, control=list(maxit=1000))

s1 <- numeric(4)
s2 <- numeric(4)
for(i in 1:4) {
  res <- optim(c(runif(1), runif(1)), matT_log_density)
  s1[i] <- res$par[1]
  s2[i] <- res$par[2]
}
cat("Sigma 1:",s1,"\n")
cat("Sigma 2:",s2,"\n")
