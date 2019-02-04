library(cluster)
library(ClusterR)
library(psych)
library(Rcpp)
source("include.R")
sourceCpp("fastCorr.cpp")

# =============================================================================
# RUN K-MEDOIDS CLUSTERING
# =============================================================================

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for pseudocount missing!\n")
}

pseudocount <- as.numeric(args[1])

if(!exists("filtered")) {
  filtered <- filter_data()
  # data is ordered by sname, then collection_date
}

subset_to <- nsamples(filtered)
K <- 500
sample_idx <- sample(nsamples(filtered))[1:subset_to]

ilr_data <- apply_ilr(filtered, pseudocount=pseudocount)
ilr_data <- ilr_data[,sample_idx] # subset
ilr_data <- t(apply(ilr_data, 1, function(x) x - mean(x)))

# PAM from ClusterR package (Rcpp)

res_cr <- Cluster_Medoids(ilr_data, K, distance_metric="pearson_correlation", verbose=FALSE, threads=4)
save(res_cr, file=paste("cluster_results_",pseudocount,".RData",sep=""))

# PAM in R

if(FALSE) {
diss_mat <- 1 - abs(fastCorr(ilr_data))
res_pam <- pam(diss_mat, K, diss=TRUE)
# alternatively, w/o dissimilarity matrix: res_pam <- pam(t(ilr_data), K)
}

if(FALSE) {
  sweep_pairs <- combn(1:length(pc_sweep), 2)
  for(pair in 1:dim(sweep_pairs)[2]) {
    pc1_idx <- sweep_pairs[1,pair]
    pc2_idx <- sweep_pairs[2,pair]
    c1 <- clusterings[[pc1_idx]]
    c2 <- clusterings[[pc2_idx]]

    pairs_compared <- 0
    pairs_matched <- 0
    for(i in 1:subset_to) {
      for(j in 1:subset_to) {
        if(i != j) {
          # for a given pair
          if(c1[i] == c1[j]) {
            # if matched on the first run
            if(c2[i] == c2[j]) {
              pairs_matched <- pairs_matched + 1
            }
            pairs_compared <- pairs_compared + 1
          }
        }
      }
    }
    cat("Combination",pc_sweep[pc1_idx],"vs.",pc_sweep[pc2_idx],":",round((pairs_matched/pairs_compared),digits=3),"\n")
  }
}
