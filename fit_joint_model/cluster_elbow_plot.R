# this file repeatedly applies clustering (PAM again) to the output of the individual model fits
#
# the idea is to get a within/between cluster error and build an elbow plot:
#   https://www.wikiwand.com/en/Elbow_method_(clustering)
# to evaluate a "good" cluster number or a set of possible "modes" associated with dynamics
#
# this error is directly spit to STDOUT, so the output of this must be piped to a file;
# that output can be parsed by `parse_PAM_elbow_scores.pl` to give {output_dir}/all_elbow_scores.txt

relative_path <- ".."

library(ClusterR)

source(file.path(relative_path,"include/R/general.R")) # for sname_list

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

# useful to use the _MAP version of these distances (WAY smaller) for testing
distmat <- readRDS(file.path(relative_path,output_dir,paste0("saved_distance_Sigma_",level,".rds")))

df <- data.frame(k=c(), ss=c())
for(k in 69:round(length(sname_list)/2)) {
  cat("Clustering into",k,"components...\n")
  sample_cr <- Cluster_Medoids(distmat, k, verbose=FALSE, threads=1)
  intra_dist <- sample_cr$silhouette_matrix$intra_clust_dissim
  # clusters with one member have intracluster distances of NaN
  intra_dist[is.na(intra_dist)] <- 0
  error <- sum(intra_dist)
  cat("\tError:",round(error,3),"\n")
  df <- rbind(df, data.frame(k=k, ss=error))
}

