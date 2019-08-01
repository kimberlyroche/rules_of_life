# plot time course and mean covariance (ILR) for a given individual

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript 04_visualize_individual.R family ACA ZIZ", call.=FALSE)
}
level <- args[1]
baboons <- c(args[2])
if(length(args) > 2) {
  for(i in 3:length(args)) {
    baboons <- c(baboons, args[i])
  }
}

source("include.R")

glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.33, collapse_level=level, verbose=TRUE)

for(baboon in baboons) {
  cat("Baboon:",baboon,", level:",level,"\n")

  indiv_data <- subset_samples(data, sname==baboon)
  cat("\tPlotting timecourse...\n")
  # these functions already prepends with 'plot' -- watch out!
  plot_timecourse_phyloseq(indiv_data, paste0("basset/",level,"/",baboon,"_timecourse"), gapped=FALSE,
                                     legend=TRUE, legend_level=level)

  # impossible to differentiate the colors here; need to truncate to top 10 probably
  #load(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.RData"))
  #collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) }))
  #propLambda <- ilrInv(collLambda) # applied with default basis, so V=NULL should be ok?
  #df <- gather_array(propLambda, proportion, sample, enzyme)
  #df$sample <- as.factor(df$sample)
  #df$enzyme <- as.factor(df$enzyme) # "enzyme" is taxon here
  #plot_timecourse_metagenomics(df, save_filename=paste0("basset/",level,"/",baboon,"_baseline"))

  cat("\tPlotting mean covariance/correlation...\n")
  # load Sigma samples
  fn <- paste0("subsetted_indiv_data/",level,"/",baboon,"_bassetfit.RData")
  load(fn)
  meanSigma <- apply(Sigma, c(1,2), mean)
  cat("\t\tTrace:",sum(diag(meanSigma)),"\n")
  df <- driver::gather_array(meanSigma, "value", "feature_row", "feature_col")
  p <- ggplot(df, aes(feature_row, feature_col)) +
         geom_tile(aes(fill = value), colour = "white") +
         scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
  ggsave(paste0("plots/basset/",level,"/",baboon,"_mean_cov.png"), plot=p, scale=1.5, width=7, height=6, units="in", dpi=72)
  meanSigma_corr <- cov2cor(meanSigma)
  df <- driver::gather_array(meanSigma_corr, "value", "feature_row", "feature_col")
  p <- ggplot(df, aes(feature_row, feature_col)) +
         geom_tile(aes(fill = value), colour = "white") +
         scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
  ggsave(paste0("plots/basset/",level,"/",baboon,"_mean_corr.png"), plot=p, scale=1.5, width=7, height=6, units="in", dpi=72)
}
