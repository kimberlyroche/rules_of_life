# this file uses PERMANOVA to ask whether distances between individuals in the space of dynamics
# explain differences in any of the fitness/stress annotations provided (i.e. born_in_drought)
# tl;dr -- it doesn't seem like any of these are convincing hits

relative_path <- ".."

library(vegan)

source(file.path(relative_path,"include/R/GP.R"))

levels <- c("phylum", "family", "genus")

for(level in levels) {
  cat("Level:",level,"\n\n")

  # below are the estimates fit by `preprocessing/formalize_parameters.R`
  # for data filtered at the given level at at least a 5-count in 20% of samples
  if(level == "phylum") {
    se_weight <- sqrt(2.052)
    per_weight <- sqrt(0.228)
    alr_ref <- 9
  }
  if(level == "family") {
    se_weight <- sqrt(2.586)
    per_weight <- sqrt(0.287)
    alr_ref <- 22
  }
  if(level == "genus") {
    se_weight <- sqrt(2.632)
    per_weight <- sqrt(0.292)
    alr_ref <- 59
  }

  baboons <- sort(sname_list)
  all_samples_Sigma <- NULL
  all_samples_Lambda <- NULL
  P <- NULL
  V <- NULL
  for(i in 1:length(baboons)) {
    baboon <- baboons[i]
    fit_obj <- fit_GP(baboon, level, se_weight=se_weight, per_weight=per_weight, wn_weight=0,
                      alr_ref=alr_ref, dd_se=90, save_append="", date_lower_limit=NULL, date_upper_limit=NULL, verbose=FALSE, mean_only=TRUE)
    dim(fit_obj$fit$Eta) <- c(nrow(fit_obj$fit$Eta), ncol(fit_obj$fit$Eta), 1)
    dim(fit_obj$fit$Lambda) <- c(nrow(fit_obj$fit$Lambda), ncol(fit_obj$fit$Lambda), 1)
    dim(fit_obj$fit$Sigma) <- c(nrow(fit_obj$fit$Sigma), ncol(fit_obj$fit$Sigma), 1)
    if(is.null(V)) {
      V <- driver::create_default_ilr_base(ncategories(fit_obj$fit))
    }
    fit.ilr <- to_ilr(fit_obj$fit, V)
    Lambda <- fit.ilr$Lambda
    Sigma <- fit.ilr$Sigma
     if(is.null(all_samples_Sigma)) {
      P <- fit.ilr$D-1
    }
  
    # note: this is calculated incorrectly in PibbleCollapse_Uncollapse.cpp
    upsilonN <- fit.ilr$upsilon + fit.ilr$N
    Sigma <- Sigma / (upsilonN - fit.ilr$D)
    Sigma <- Sigma / (upsilonN - P - 1)
    # this should be the correct invWishart mean
    if(is.null(all_samples_Sigma)) {
      all_samples_Sigma <- matrix(NA, P, P*length(baboons))
    }
    all_samples_Sigma[,((i-1)*P+1):(i*P)] <- Sigma
    if(is.null(all_samples_Lambda)) {
      all_samples_Lambda <- matrix(NA, length(baboons), P)
    }
    collLambda <- apply(Lambda[,,1], 1, mean) # P x 1
    all_samples_Lambda[i,] <- collLambda
  }
  distance_mat <- Riemann_dist_samples(all_samples_Sigma, length(baboons), 1)

  df <- readRDS(file.path(relative_path,plot_dir,"basset",paste0(level,"/Sigma_ordination_centroids.rds")))
  glom_data <- load_glommed_data(level=level, replicates=TRUE)

  x_factors <- list(c("group", 1),
            c("matgroup", 1), 
            c("mom", 1),
            c("dad", 1), 
            c("momrank", 1), 
            c("drought", 1), 
            c("largegroup", 1), 
            c("momdied", 1), 
            c("competingsib", 1), 
            c("earlyadversity", 1), 
            c("birthrate_all", 1), 
            c("birthrate_surviving", 1),
            c("counts", 0),
            c("density", 0))

  for(x_factor in x_factors) {
    cat(x_factor[1],"\t")
    labeled_df <- get_other_labels(df, glom_data, unique(df$labels), annotation=x_factor[1])
    labels <- labeled_df[rownames(labeled_df) %in% baboons,]$labels
    sub_distance_mat <- distance_mat[!is.na(labels),!is.na(labels)]
    labels <- labels[!is.na(labels)]
    labels <- data.frame(label=labels)
    obj <- adonis(sub_distance_mat ~ label, data=labels, permutations=10000)
    R2 <- obj$aov.tab$R2
    Pval <- obj$aov.tab$`Pr(>F)`
    cat(round(R2[1], 3),"\t",round(Pval[1], 3),"\n")
  }
}
