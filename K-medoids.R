library(cluster)
library(ClusterR)
library(psych)
library(Rcpp)
source("include.R")

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

#subset_to <- nsamples(filtered)
subset_to <- 200
#K <- 500
K <- 5
sample_idx <- sample(nsamples(filtered))[1:subset_to]

# use PAM from ClusterR package (Rcpp)
if(TRUE) {
  counts <- otu_table(filtered)@.Data + pseudocount
  # dissimilarity as Aitchison distance
  diss_mat <- matrix(0, nrow=subset_to, ncol=subset_to)
  for(i in 1:subset_to) {
    for(j in 1:subset_to) {
      if(j >= i) {
        x.i <- as.vector(counts[i,])
        x.j <- as.vector(counts[j,])
        diss_mat[i,j] <- aitchison_dist(x.i, x.j)
      } else {
        diss_mat[i,j] <- diss_mat[j,i]
      }
    }
  }
  res_cr <- Cluster_Medoids(diss_mat, K, verbose=FALSE, threads=4)
} else {
  ilr_data <- apply_ilr(filtered, pseudocount=pseudocount)
  ilr_data <- ilr_data[,sample_idx] # subset
  ilr_data <- apply(ilr_data, 1, function(x) x - mean(x))
  res_cr <- Cluster_Medoids(ilr_data, K, distance_metric="pearson_correlation", verbose=FALSE, threads=4)
}
save(res_cr, file=paste("testing_cluster_results_",pseudocount,".RData",sep=""))

# PAM in R

if(FALSE) {
ilr_data <- t(ilr_data)
diss_mat <- 1 - abs(fastCorr(ilr_data))
res_pam <- pam(diss_mat, K, diss=TRUE)
# alternatively, w/o dissimilarity matrix: res_pam <- pam(t(ilr_data), K)
}
