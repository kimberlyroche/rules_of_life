library(cluster)
library(ClusterR)
library(psych)
library(Rcpp)

source("include.R")

# =============================================================================
# run PAM clustering (K-medoids)
# =============================================================================

# usage: pipeline_PAM.R 0.65

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for pseudocount missing!\n")
}

pseudocount <- as.numeric(args[1])

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)

# subset to so many samples and partition into K clusters
subset_to <- 200 # testing
k <- 5
sample_idx <- sample(nsamples(filtered))[1:subset_to]

# use PAM from ClusterR package (Rcpp)
counts <- otu_table(filtered)@.Data + pseudocount
diss_mat <- matrix(0, nrow=subset_to, ncol=subset_to)
for(i in 1:subset_to) {
  for(j in 1:subset_to) {
    if(j >= i) {
      x.i <- as.vector(counts[i,])
      x.j <- as.vector(counts[j,])
      # code.base::dist() also includes Aitchison distance
      diss_mat[i,j] <- aitchison_dist(x.i, x.j)
    } else {
      diss_mat[i,j] <- diss_mat[j,i]
    }
  }
}
res_cr <- Cluster_Medoids(diss_mat, k, verbose=FALSE, threads=4)
save(res_cr, file=paste("PAM_results_pc",pseudocount,".RData",sep=""))
