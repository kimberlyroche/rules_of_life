#library(stray)
devtools::load_all("/data/mukherjeelab/labraduck")

library(phyloseq)
library(dplyr)
library(driver)
library(tidyverse)

# X is Q x N as in other kernels
# bandwidth: rho as chosen gives antiphase observations a correlation of ~0.1
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

get_predictions <- function(X, fit, n_samples=2000){
  cat("Predicting from 1 to",max(X),"\n")
  X_predict <- t(1:(max(X))) # time point range, fill in any missing
  predicted <- predict(fit, X_predict, iter=n_samples) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

fit_to_baboon <- function(baboon, indiv_data, Gamma, date_lower_limit=NULL, date_upper_limit=NULL) {
  Y_full <- indiv_data$ys
  observations_full <- indiv_data$observation_vec
  
  if(!is.null(date_lower_limit) & !is.null(date_upper_limit)) {
    # require both be present for now
    min_idx <- min(which(names(observations_full) >= date_lower_limit))
    max_idx <- max(which(names(observations_full) <= date_upper_limit))
    
    # subset dates
    Y <- t(Y_full[min_idx:max_idx,])
    colnames(Y) <- NULL; rownames(Y) <- NULL
    # set first observation to t=1
    observations <- matrix(observations_full[min_idx:max_idx], nrow=1) - observations_full[min_idx] + 1
  } else {
    Y <- t(Y_full)
    colnames(Y) <- NULL
    rownames(Y) <- NULL
    observations <- matrix(observations_full, nrow=1)
  }

  D <- nrow(Y)
  N <- ncol(Y)
  
  # ALR prior covariance
  upsilon <- D-1+10 # lesser certainty
  # supsilon <- D-1+20 # greater certainty; this should tighten the distribution around this mean
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference;
  # take diag as covariance over log abundances
  Xi <- GG%*%(diag(D)*1)%*%t(GG)
  # mean-center
  Xi <- Xi*(upsilon-D-1)
  
#  alr_ys <- driver::alr(t(Y) + 0.5)
#  alr_means <- colMeans(alr_ys)
#  Theta <- function(X) matrix(alr_means, D-1, ncol(X))

  Theta <- function(X) matrix(10, D-1, ncol(X))

  fit <- stray::basset(Y, observations, upsilon, Theta, Gamma, Xi)
  #fit.clr <- to_clr(fit)
  fit.alr <- to_alr(fit, D)
  return(list(Y=Y, X=observations, fit=fit.alr))
}

plot_predictions <- function(fit_obj, predict_obj, LR_coord=1, save_name="out") {
  Y.alr <- driver::alr(t(fit_obj$Y) + 0.5)
  observations <- fit_obj$X
  alr_tidy <- gather_array(Y.alr, "LR_value", "timepoint", "LR_coord")
  
  # replace timepoints with observation dates
  map <- data.frame(timepoint=1:length(observations), observation=c(observations))
  alr_tidy <- merge(alr_tidy, map, by="timepoint")
  alr_tidy <- alr_tidy[,!(names(alr_tidy) %in% c("timepoint"))]

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
    geom_point(data=alr_tidy[alr_tidy$LR_coord==LR_coord,], aes(x=observation, y=LR_value), alpha=0.5) +
    theme_minimal() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle=45))
  if(LR_coord == 1) {
    p <- p + ylab("ALR(Bifidobacteriaceae/Helicobacteraceae)")
  } else if(LR_coord == 2) {
    p <- p + ylab("ALR(class Kiritimatiellae/Helicobacteraceae)")
  } else if(LR_coord == 7) {
    p <- p + ylab("ALR(Muribaculaceae/Helicobacteraceae)")
  } else if(LR_coord == 19) {
    p <- p + ylab("ALR(order Mollicutes RF39/Helicobacteraceae)")
  } else {
    p <- p + ylab("ALR coord")
  }
  ggsave(paste0("plots/basset/",save_name,".png"), scale=2, width=12, height=3, units="in", dpi=100)
}

# for reference, individuals passable as arguments are:
# "DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI"

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Usage: Rscript basset_walkthrough.r ACA", call.=FALSE)
}

baboon <- args[1]

load(paste0("subsetted_indiv_data/",baboon,"_data.RData"))

#Gamma <- function(X) PER(X, period=365) # periodic only
Gamma <- function(X) SE(X, sigma=1, rho=100, jitter=1e-8) # squared exponential only
#Gamma <- function(X) 0.15*PER(X, sigma=1, rho=20, period=365, jitter=0) + 0.85*SE(X, sigma=1, rho=100, jitter=0) + (1e-8)*diag(ncol(X)) # pretty arbitrary
# will additive combinations like this always be PSD?

#fit_obj <- fit_to_baboon(baboon, indiv_data, Gamma, date_lower_limit="2001-10-01", date_upper_limit="2003-11-30")
fit_obj <- fit_to_baboon(baboon, indiv_data, Gamma)
predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=1000) # interpolates
plot_predictions(fit_obj, predict_obj, LR_coord=1, save_name=paste0(baboon,"_LR1_2"))
plot_predictions(fit_obj, predict_obj, LR_coord=2, save_name=paste0(baboon,"_LR2_2"))
plot_predictions(fit_obj, predict_obj, LR_coord=7, save_name=paste0(baboon,"_LR7_2"))
plot_predictions(fit_obj, predict_obj, LR_coord=19, save_name=paste0(baboon,"_LR19_2"))

#save(fit_obj, file=paste0("bassetfit_",baboon,".RData"))





