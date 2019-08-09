library(stray)
library(driver)

source("include.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Testing usage: Rscript 05_visualize_extrema.R family", call.=FALSE)
}
level <- args[1]

plot_extreme_Lambda <- function(coordinate, no_indiv=10, trundate_taxa=10, save_filename="../test") {
  df_centroids <- readRDS("plots/basset/family/Lambda_ordination_centroids.rds")
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))

  # first pass: get non-tiny proportions to retain
  allLambda <- NULL
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"(1)\n")
    fit_obj <- readRDS(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.rds"))
    fit.clr <- to_clr(fit_obj$fit)
    Lambda <- fit.clr$Lambda
    collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) }))
    propLambda <- clrInv(collLambda)
    if(is.null(allLambda)) {
      allLambda <- propLambda
    } else {
      allLambda <- rbind(allLambda, propLambda)
    }
  }
  mean_all <- apply(allLambda, 2, mean)
  retain_prop <- mean_all > 0.01
  for(i in 1:length(retain_prop)) {
    if(retain_prop[i]) {
      cat("Retaining CLR coord",i,"with mean",mean_all[i],"\n")
    }
  }

  # stupid feature naming is for the benefit of plot_timecourse_metagenomics
  df <- data.frame(sample=c(), enzyme=c(), proportion=c())

  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"(2)\n")
    fit_obj <- readRDS(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.rds"))
    fit.clr <- to_clr(fit_obj$fit)
    Lambda <- fit.clr$Lambda
    collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) }))
    #propLambda <- ilrInv(collLambda, V=V) # applied with default basis, so V=NULL should be ok?
    propLambda <- clrInv(collLambda)
    avgProp <- colMeans(propLambda)[retain_prop] # truncate to make readable
    df <- rbind(df, data.frame(sample=rep(baboon, length(avgProp)),
                               enzyme=as.factor(1:length(avgProp)), proportion=avgProp))
  }
  df$sample <- as.factor(df$sample)
  plot_timecourse_metagenomics(df, save_filename=paste0(save_filename))
}

plot_extreme_Sigma <- function(coordinate, no_indiv=10, save_filename="test") {
  df_centroids <- readRDS("plots/basset/family/Sigma_ordination_centroids.rds")
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))

  df_on <- data.frame(sample=c(), feature=c(), value=c())
  df_off <- data.frame(sample=c(), feature=c(), value=c())
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"\n")
    Sigma <- readRDS(paste0("subsetted_indiv_data/family/",baboon,"_bassetfit.rds"))$fit$Sigma
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

plot_diag_Sigma <- function(lrtransform="alr", save_filename="test") {
  indiv_obj <- fitted_individuals(level="family")
  individuals <- indiv_obj$individuals
  pattern_str <- indiv_obj$pattern_str
  regexpr_str <- indiv_obj$regexpr_str

  fit_obj <- readRDS(paste0("subsetted_indiv_data/",level,"/",individuals[1],regexpr_str))
  D <- nrow(fit_obj$Y)
  P <- D
  if(lrtransform == "alr" | lrtransform == "ilr") {
    P <- D-1
  }

  msd_mat <- matrix(0, length(individuals), P)
  tax_var_mat <- matrix(0, length(individuals), P)
  for(i in 1:length(individuals)) {
    baboon <- individuals[i]
    cat("Loading fit for",baboon,"\n")
    fit <- readRDS(paste0("subsetted_indiv_data/",level,"/",baboon,regexpr_str))$fit
    clr_ys <- driver::clr(t(fit$Y) + 0.5)
    tax_var <- apply(clr_ys, 2, sd)
    if(lrtransform == "clr") {
      fit <- to_clr(fit)
    } else if(lrtransform == "ilr") {
      V <- driver::create_default_ilr_base(ncategories(fit))
      fit <- to_ilr(fit, V)
    }
    Sigma <- fit$Sigma
    meanSigmaDiag <- diag(apply(Sigma, c(1,2), mean))
    msd_mat[i,] <- meanSigmaDiag
    tax_var_mat[i,] <- tax_var
  }
  indiv_dist <- dist(msd_mat)
  clust_order <- hclust(indiv_dist)$order
  msd_ordered <- msd_mat[clust_order,]
  tax_var_mean <- apply(tax_var_mat, 2, mean)
  df <- gather_array(msd_ordered)

  individuals <- individuals[clust_order]

  df <- as.data.frame(t(msd_ordered))
  colnames(df) <- individuals
  df_long <- gather_array(df, "value", "taxon", "baboon")
  df_long$baboon <- as.factor(df_long$baboon)
  levels(df_long$baboon) <- individuals

  p <- ggplot(df_long, aes(baboon, taxon)) +
         geom_tile(aes(fill = value), colour = "white") +
         scale_fill_gradient2(low = "white", high = "darkred") +
         theme(axis.text.x = element_text(face="bold", size=17, angle=90))
  ggsave(paste0(save_filename,"_",lrtransform,".png"), plot=p, scale=2, width=16, height=8, units="in", dpi=72)
  cat("Mean empirical LR variance:\n")
  for(p in 1:length(tax_var_mean)) {
    cat("\t",tax_var_mean[p],"\n")
  }
}

plot_extreme_Lambda(coordinate="mean_x", no_indiv=10, trundate_taxa=10, save_filename=paste0("basset/",level,"/Lambda_ordination_PC1_extrema"))
plot_extreme_Lambda(coordinate="mean_y", no_indiv=10, trundate_taxa=10, save_filename=paste0("basset/",level,"/Lambda_ordination_PC2_extrema"))

#plot_extreme_Sigma("mean_x", no_indiv=10, save_filename=paste0("plots/basset/",level,"/Sigma_ordination_PC1_extrema"))
#plot_extreme_Sigma("mean_y", no_indiv=10, save_filename=paste0("plots/basset/",level,"/Sigma_ordination_PC2_extrema"))

#plot_diag_Sigma(lrtransform="alr", save_filename=paste0("plots/basset/",level,"/Sigma_diagonal"))
#plot_diag_Sigma(lrtransform="clr", save_filename=paste0("plots/basset/",level,"/Sigma_diagonal"))
#plot_diag_Sigma(lrtransform="ilr", save_filename=paste0("plots/basset/",level,"/Sigma_diagonal"))

# transform through all possible ALR references
# look for low-variance ratios across individuals

# to ILR
# V <- driver::create_default_ilr_base(ncategories(fit_obj$fit))
# fit.ilr <- to_ilr(fit_obj$fit, V)

# to CLR
# fit.clr <- to_clr(fit_obj$fit)
