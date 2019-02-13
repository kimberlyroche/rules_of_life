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
ilr_data <- apply(ilr_data, 1, function(x) x - mean(x))

# PAM from ClusterR package (Rcpp)

res_cr <- Cluster_Medoids(ilr_data, K, distance_metric="pearson_correlation", verbose=FALSE, threads=4)
save(res_cr, file=paste("cluster_results_",pseudocount,".RData",sep=""))

# PAM in R

if(FALSE) {
ilr_data <- t(ilr_data)
diss_mat <- 1 - abs(fastCorr(ilr_data))
res_pam <- pam(diss_mat, K, diss=TRUE)
# alternatively, w/o dissimilarity matrix: res_pam <- pam(t(ilr_data), K)
}
