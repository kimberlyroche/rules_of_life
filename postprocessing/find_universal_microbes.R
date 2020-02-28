# this file represents one attempt to identify "universal" microbes (those will a set of associations to other
# taxa that is preserved [the set of relationships is preserved, that is] across hosts)
#
# note: need to re-run this and decide whether this approach is worth keeping

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

level <- "family"

if(FALSE) {
  
  # pull all fitted models
  models <- get_fitted_modellist_details(level=level, MAP=TRUE)
  
  Sigmas <- list()
  
  for(i in 1:length(models$individuals)) {
    host <- models$individuals[i]
    cat("Processing",host,"\n")
    fit <- readRDS(models$model_list[i])$fit
    # fit "incorrect number of dimensions error"
    dim(fit$Eta) <- c(dim(fit$Eta), 1)
    dim(fit$Lambda) <- c(dim(fit$Lambda), 1)
    dim(fit$Sigma) <- c(dim(fit$Sigma), 1)
    fit.clr <- to_clr(fit)
    Sigmas[[host]] <- fit.clr$Sigma[,,1]
  }
  
  # build some readable labels
  ref <- readRDS(file.path(relative_path,data_dir,"filtered_family_5_20.rds"))
  tax <- tax_table(ref)@.Data
  labels <- as.vector(tax[,5])
  for(i in 1:length(labels)) {
    if(i == 42) {
      # this is "Other" category for family-level data filtered at 5-counts in 20% or more samples
      labels[i] <- "CLR(Other)"
    } else {
      if(is.na(labels[i])) {
        for(j in 4:1) {
          if(!is.na(tax[i,j])) {
            if(j == 1) {
              labels[i] <- paste0("CLR(kingdom ",tax[i,j],")")
            }
            if(j == 2) {
              labels[i] <- paste0("CLR(phylum ",tax[i,j],")")
            }
            if(j == 3) {
              labels[i] <- paste0("CLR(class ",tax[i,j],")")
            }
            if(j == 4) {
              labels[i] <- paste0("CLR(order ",tax[i,j],")")
            }
            break
          }
        }
      } else {
        labels[i] <- paste0("CLR(family ",labels[i],")")
      }
    }
  }
  
  # build a map from which we can translate an interaction number with the original pair of taxa
  # (I'm sure there's a more elegant way to do this)
  label_pairs <- matrix(NA, models$D, models$D)
  for(i in 1:models$D) {
    for(j in 1:models$D) {
      if(i < j) {
        label_pairs[i,j] <- paste0(i,"_",j)
      }
    }
  }
  label_pairs <- label_pairs[upper.tri(label_pairs, diag=F)]
  
  # plot all in heatmap as hosts x pairs of microbial interaction
  interaction_no <- (models$D^2)/2 - models$D/2
  interaction_pairs <- matrix(NA, length(models$individuals), interaction_no)
  for(i in 1:length(models$individuals)) {
    host <- models$individuals[i]
    # convert this host's MAP covariance to correlation
    corr_mat <- cov2cor(Sigmas[[host]])
    # stack all unique pairwise correlations between microbes in a row associated with this host
    interaction_pairs[i,] <- corr_mat[upper.tri(Sigmas[[host]], diag=F)]
  }
  
  # get distance of each pair
  d <- dist(t(interaction_pairs))
  clustering <- hclust(d)
  
  # reorder
  interaction_pairs <- interaction_pairs[,clustering$order]
  
  avg_interaction <- colMeans(interaction_pairs)
  
  # get strong-ish average interactions
  interesting_idx <- which(abs(avg_interaction) >= 0.4)
  for(p_idx in interesting_idx) {
    interesting_pair <- label_pairs[p_idx]
    microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
    cat("Interesting pair (",p_idx,"):",labels[microbe_pair[1]],",",labels[microbe_pair[2]],"\n")
  }
  
  # plot heatmap
  df <- gather_array(interaction_pairs, "correlation", "host", "pair")
  p <- ggplot(df, aes(pair, host)) +
    geom_tile(aes(fill = correlation), colour = "white") +
    scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave(file.path(relative_path,plot_dir,"microbe_pair_correlations.png"), p, units="in", dpi=150, height=5, width=15)
  
  # sample across some random hosts; do these taxa actually appear correlated at the level of log relative abundances?
  microbe_pair <- as.numeric(strsplit(label_pairs[1], "_")[[1]])
  sampled_hosts <- models$individuals[sample(1:length(models$individuals))[1:2]]
  # fits <- list()
  # for(host in sampled_hosts) {
  #   # get full posteriors
  #   fits[[host]] <- readRDS(file.path(relative_path,model_dir,level,paste0(host,"_bassetfit.rds")))
  #   # sanity check universal pairs via predictive plots
  #   plot_ribbons_individuals(c("CRU"), level, timecourse=FALSE, covcor=FALSE, predict_coords=c(1))
  # }
  
}


