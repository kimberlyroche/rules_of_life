library(cluster)
library(ClusterR)
library(psych)
library(Rcpp)
source("include.R")
sourceCpp("fastCorr.cpp")

# =============================================================================
# RUN K-MEDOIDS CLUSTERING
# =============================================================================

pseudocount <- 0.1
cat("Reading in counts\n")
load("glom_data_genus.RData")
filtered <- filter_data(data=glom_data)
counts <- t(otu_table(filtered)) + pseudocount
rm(filtered)

testing <- TRUE

if(testing) {
  cat("Loading clustering results for PC=0.1\n")
  load("testing_cluster_results_0.1.RData")
  res_1 <- res_cr
  cat("Loading clustering results for PC=0.65\n")
  load("testing_cluster_results_0.65.RData")
  res_2 <- res_cr
  cat("Loading clustering results for PC=1.0\n")
  load("testing_cluster_results_1.RData")
  res_3 <- res_cr
  rm(res_cr)
} else {
  cat("Loading clustering results for PC=0.1\n")
  load("cluster_results_0.1.RData")
  res_1 <- res_cr
  cat("Loading clustering results for PC=0.65\n")
  load("cluster_results_0.65.RData")
  res_2 <- res_cr
  cat("Loading clustering results for PC=1.0\n")
  load("cluster_results_1.RData")
  res_3 <- res_cr
  rm(res_cr)
}

labels <- c("0.1", "0.65", "1.0")
clusterings <- list(res_1, res_2, res_3)
sweep_pairs <- combn(1:length(clusterings), 2)
subset_to <- sample(length(clusterings[[1]]$clusters))[1:200]

for(pair in 1:dim(sweep_pairs)[2]) {
  pc1_idx <- sweep_pairs[1,pair]
  pc2_idx <- sweep_pairs[2,pair]
  c1 <- clusterings[[pc1_idx]]$clusters
  c2 <- clusterings[[pc2_idx]]$clusters

  pairs_compared <- 0
  pairs_matched <- 0
  for(i in subset_to) {
    for(j in subset_to) {
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
  cat("Combination",labels[pc1_idx],"vs.",labels[pc2_idx],":",round((pairs_matched/pairs_compared),digits=5),"\n")
}

within <- c()
between <- c()
assignments <- res_1$clusters[subset_to]
for(i in 1:(length(assignments)-1)) {
  for(j in (i+1):length(assignments)) {
    c.i <- as.vector(counts[,i])
    c.j <- as.vector(counts[,j])
    d <- aitchison_dist(c.i, c.j)
    if(assignments[i] == assignments[j]) {
      within[length(within)+1] <- d
    } else {
      between[length(between)+1] <- d
    }
  }
}

png("hist_within.png")
hist(within)
dev.off()

png("hist_between.png")
hist(between)
dev.off()

cat("Average (sampled) within-cluster distance:",mean(within),"\n")
cat("Average (sampled) between-cluster distance:",mean(between),"\n")

