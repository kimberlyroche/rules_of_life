source("include.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Testing usage: Rscript 05_visualize_vectors.R family", call.=FALSE)
}
level <- args[1]

plot_extreme_Lambda <- function(coordinate, no_indiv=10, trundate_taxa=10, save_filename="../test") {
  load("plots/basset/family/Lambda_ordination_centroids.RData")
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))

  # stupid feature naming is for the benefit of plot_timecourse_metagenomics
  df <- data.frame(sample=c(), enzyme=c(), proportion=c())

  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"\n")
    load(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.RData"))
    collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) }))
    propLambda <- ilrInv(collLambda, V=V) # applied with default basis, so V=NULL should be ok?
    avgProp <- colMeans(propLambda)[1:truncate_taxa] # truncate to make readable
    df <- rbind(df, data.frame(sample=rep(baboon, truncate_taxa), enzyme=as.factor(1:truncate_taxa), proportion=avgProp))
  }
  df$sample <- as.factor(df$sample)
  plot_timecourse_metagenomics(df, save_filename=paste0(save_filename))
}

plot_extreme_Sigma <- function(coordinate, no_indiv=10, save_filename="../test") {
  load("plots/basset/family/Sigma_ordination_centroids.RData")
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))

  df_on <- data.frame(sample=c(), feature=c(), value=c())
  df_off <- data.frame(sample=c(), feature=c(), value=c())
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"\n")
    load(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.RData"))
    meanSigma <- apply(Sigma, c(1,2), mean)
    diagSigma <- diag(meanSigma)
    df_on <- rbind(df_on, data.frame(sample=rep(baboon, length(diagSigma)), feature=as.factor(1:length(diagSigma)), value=diagSigma))
    upperSigma <- c(meanSigma[upper.tri(meanSigma)])
    df_off <- rbind(df_off, data.frame(sample=rep(baboon, length(upperSigma)), feature=as.factor(1:length(upperSigma)), value=upperSigma))
  }
  df_on$sample <- as.factor(df_on$sample)
  df_off$sample <- as.factor(df_off$sample)

  p <- ggplot(df_on, aes(feature, sample)) +
         geom_tile(aes(fill = value), colour = "white") +
         scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  ggsave(paste0(save_filename,"_ondiag.png"), plot=p, scale=2, width=6, height=6, units="in", dpi=72)
  p <- ggplot(df_off, aes(feature, sample)) +
         geom_tile(aes(fill = value), colour = "white") +
         scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  ggsave(paste0(save_filename,"_offdiag.png"), plot=p, scale=2, width=10, height=6, units="in", dpi=72)
}

#plot_extreme_Lambda("mean_x", no_indiv=10, trundate_taxa=10,
#                    save_filename=paste0("basset/",level,"/Lambda_ordination_PC1_extrema"))
#plot_extreme_Lambda("mean_y", no_indiv=10, trundate_taxa=10,
#                    save_filename=paste0("basset/",level,"/Lambda_ordination_PC2_extrema"))

#plot_extreme_Sigma("mean_x", no_indiv=10, save_filename=paste0("plots/basset/",level,"/Sigma_ordination_PC1_extrema"))
plot_extreme_Sigma("mean_y", no_indiv=10, save_filename=paste0("plots/basset/",level,"/Sigma_ordination_PC2_extrema"))


