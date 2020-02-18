library(ClusterR)
source("include/R/general.R") # for sname_list

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

#distmat <- readRDS(paste0(output_dir,"saved_distance_Sigma_",level,"_MAP.rds"))
distmat <- readRDS(paste0(output_dir,"saved_distance_Sigma_",level,".rds"))

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

quit()

p <- ggplot(df) +
      geom_point(aes(x=k, y=ss), size=2) +
      ylim(0, max(df$ss)) +
      theme(legend.position = "none")
ggsave(paste0("elbow_intracluster_",level,".png"), units="in", dpi=100, height=10, width=10)

# visualize
embedding <- cmdscale(distmat)
df <- data.frame(x=embedding[,1], y=embedding[,2], sname=sname_list)
p <- ggplot(df) +
      geom_point(aes(x=x, y=y, color=sname), size=2) +
      theme(legend.position = "none")
ggsave(paste0("elbow_indiv_",level,".png"), units="in", dpi=100, height=10, width=10)

