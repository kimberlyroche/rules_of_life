library(cluster)
library(psych)
library(Rcpp)
source("include.R")
sourceCpp("fastCorr.cpp")

pdf(NULL)

# =============================================================================
# GET REPRODUCIBILITY SCORE BETWEEN CLUSTERING RUNS
# =============================================================================

match_runs <- function(a1, a2) {
  if(length(a1) != length(a2)) {
    return(Inf)
  }
  K <- length(unique(a1)) # works as long as clusters are assigned 1:K
  S <- matrix(0, nrow=K, ncol=K)
  # how to vectorize?
  for(k in 1:K) {
    row_occur <- a2[which(a1==k)]
    incr <- as.vector(seq(1,K))
    incr <- as.vector(lapply(incr, function(x) { occur <- length(which(row_occur==x)); if(occur > 0) { occur } else { 0 } } ))
    S[k,] <- unlist(incr)
  }
  return(S)
}

permute_cols <- function(p, a) {
  a1 <- a[[1]]
  a2 <- a[[2]]
  S <- match_runs(a1, a2)
  print(p)
  print(S)
  S <- S[,p]
  return(-tr(S))
}

# WOULD NEED TO FIND A WAY TO OPTIMIZE OVER AN ORDERING OF COLUMNS
# WHAT'S AN ALTERNATIVE WAY WE CAN JUDGE HOW WELL A CLUSTERING REPLICATES?

#a1 <- as.list(res$clustering)
#a2 <- as.list(res2$clustering)
#K <- length(unique(a1))
#res <- optim(par=as.vector(seq(1, K)), fn=permute_cols, a=list(a1, a2))

# =============================================================================
# RUN K-MEDOIDS CLUSTERING
# =============================================================================

if(!exists("filtered")) {
  filtered <- filter_data()
  # data is ordered by sname, then collection_date
}

subset_to <- 100
K <- 10
diss <- TRUE

# ILR-transform
ilr_data <- apply_ilr(filtered)
if(subset_to > 0) {
  ilr_data <- ilr_data[,sample(dim(ilr_data)[2])[1:subset_to]]
}
ilr_data <- t(apply(ilr_data, 1, function(x) x - mean(x)))

if(diss) {
  diss_mat <- 1 - abs(fastCorr(t(ilr_data)))
}

if(diss) {
  res <- pam(diss_mat, K, diss=TRUE)
} else {
  res <- pam(t(ilr_data), K)
}
cat("Results (run 1):\n")
c1 <- as.vector(res$clustering)
print(c1)

if(diss) {
  res2 <- pam(diss_mat, K, diss=TRUE)
} else {
  shuffled_idx <- sample(dim(subset_sample)[2])
  res2 <- pam(t(ilr_data[,shuffled_idx]), K)
}
cat("Results (run 2):\n")
c2 <- as.vector(res2$clustering)
print(c2)

S <- match_runs(c1, c2)
print(S)
