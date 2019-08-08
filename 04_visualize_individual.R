# plot time course and mean covariance (ILR) for a given individual

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) {
  stop("Usage: Rscript 04_visualize_individual.R family TRUE ACA ZIZ", call.=FALSE)
}
level <- args[1]
do_predict <- as.logical(args[2])
baboons <- c(args[3])
if(length(args) > 3) {
  for(i in 4:length(args)) {
    baboons <- c(baboons, args[i])
  }
}

library(stray)
source("include.R")

get_predictions <- function(X, fit, n_samples=100) {
  cat("Predicting from 1 to",max(X),"\n")
  X_predict <- t(1:(max(X))) # time point range, fill in any missing
  fit.clr <- to_clr(fit)
  predicted <- predict(fit.clr, X_predict, iter=n_samples) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

plot_predictions <- function(fit_obj, predict_obj, LR_coord=1, save_name=NULL) {
  observations <- fit_obj$X
  clr_ys <- driver::clr(t(fit_obj$Y) + 0.5)
  lr_tidy <- gather_array(clr_ys, "LR_value", "timepoint", "LR_coord")

  # replace timepoints with observation dates
  map <- data.frame(timepoint=1:length(observations), observation=c(observations))
  lr_tidy <- merge(lr_tidy, map, by="timepoint")
  lr_tidy <- lr_tidy[,!(names(lr_tidy) %in% c("timepoint"))]

  no_samples <- dim(predict_obj$Y_predict)[3]
  posterior_samples <- gather_array(predict_obj$Y_predict[LR_coord,,], "LR_value", "observation", "sample_no")
  # get quantiles

  post_quantiles <- posterior_samples %>%
    group_by(observation) %>%
    summarise(p2.5 = quantile(LR_value, prob=0.025),
              p5 = quantile(LR_value, prob=0.05),
              p10 = quantile(LR_value, prob=0.1),
              p25 = quantile(LR_value, prob=0.25),
              p50 = quantile(LR_value, prob=0.5),
              mean = mean(LR_value),
              p75 = quantile(LR_value, prob=0.75),
              p90 = quantile(LR_value, prob=0.9),
              p95 = quantile(LR_value, prob=0.95),
              p97.5 = quantile(LR_value, prob=0.975)) %>%
    ungroup()

  p <- ggplot(post_quantiles, aes(x=observation, y=mean)) +
    geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(color="blue") +
    geom_point(data=lr_tidy[lr_tidy$LR_coord==LR_coord,], aes(x=observation, y=LR_value), alpha=0.5) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45)) +
    ylab("LR coord") +
    ylim(c(-10, 10))
  if(is.null(save_name)) {
    show(p)
  } else {
    ggsave(paste0("plots/basset/",level,"/",save_name,".png"), scale=2, width=12, height=2, units="in", dpi=100)
  }
}

glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.33, collapse_level=level, verbose=TRUE)

for(baboon in baboons) {
  cat("Baboon:",baboon,", level:",level,"\n")

  indiv_data <- subset_samples(data, sname==baboon)
#  cat("\tPlotting timecourse...\n")
  # these functions already prepends with 'plot' -- watch out!
#  plot_timecourse_phyloseq(indiv_data, paste0("basset/",level,"/",baboon,"_timecourse"), gapped=FALSE,
#                                     legend=TRUE, legend_level=level)

  cat("\tPlotting mean covariance/correlation...\n")
  # load Sigma samples
  fit_obj <- readRDS(paste0("subsetted_indiv_data/",level,"/",baboon,"_bassetfit.rds"))
  Sigma <- fit_obj$fit$Sigma
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

  if(do_predict) {
    cat("\tGenerating predictions...\n")

    # a dumb hack; these need to be global
    se_weight <<- fit_obj$kernelparams$se_weight
    se_sigma <<- fit_obj$kernelparams$se_sigma
    rho_se <<- fit_obj$kernelparams$rho_se
    per_weight <<- fit_obj$kernelparams$per_weight
    per_sigma <<- fit_obj$kernelparams$per_sigma
    rho_per <<- fit_obj$kernelparams$rho_per
    wn_weight <<- fit_obj$kernelparams$wn_weight
    period <<- fit_obj$kernelparams$period

    predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=100)

    LR_coords <- c(2, 4, 27)
    #if(level == "phylum") {
    #  LR_coords <- c(1,2,3,7,8)
    #} else if(level == "family") {
    #  LR_coords <- c(1,2,3,20,21)
    #}
    for(coord in LR_coords) {
      plot_predictions(fit_obj, predict_obj, LR_coord=coord, save_name=paste0(baboon,"_",coord))
    }
  }
}
