#library(stray)
devtools::load_all("/data/mukherjeelab/labraduck")

# for reference, individuals passable as arguments are:
# "DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI"

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript 01_fit_basset.R ACA family 9", call.=FALSE)
}
baboon <- args[1]
level <- args[2]
if(length(args) > 2) {
  alr_ref <- as.numeric(args[3])
} else {
  alr_ref <- NULL
}

# testing
#date_lower_limit <- "2001-10-01"
#date_upper_limit <- "2003-11-30"
date_lower_limit <- NULL
date_upper_limit <- NULL

vizualization <- FALSE

library(phyloseq)
library(dplyr)
library(driver)
library(tidyverse)

source("include.R")

# X is Q x N as in other kernels
# bandwidth: rho as chosen gives antiphase observations a correlation of ~0.1
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

# constant kernel; has the predictable effect of flattening the whole trajectory
CONST <- function(X, sigma=1) {
  return(sigma^2)
}

get_predictions <- function(X, fit, n_samples=2000){
  cat("Predicting from 1 to",max(X),"\n")
  X_predict <- t(1:(max(X))) # time point range, fill in any missing
  predicted <- predict(fit, X_predict, iter=n_samples) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

fit_to_baboon <- function(baboon, Y, observations, Gamma, alr_ref=NULL) {
  D <- nrow(Y)
  N <- ncol(Y)
  
  # stray uses the D^th element as the ALR reference by default
  # do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }
  
  # ALR prior covariance
  upsilon <- D-1+10 # lesser certainty
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference;
  # take diag as covariance over log abundances
  Xi <- GG%*%(diag(D))%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
  
  alr_ys <- driver::alr((t(Y)+0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  #Theta <- function(X) matrix(0, D-1, ncol(X))
  
  fit <- stray::basset(Y, observations, upsilon, Theta, Gamma, Xi)
  #fit.clr <- to_clr(fit)
  #fit.alr <- to_alr(fit, alr_ref)
  return(list(Y=Y, alr_ys=alr_ys, X=observations, fit=fit))
}

plot_predictions <- function(fit_obj, predict_obj, LR_coord=1, save_name="out") {
  observations <- fit_obj$X
  alr_tidy <- gather_array(fit_obj$alr_ys, "LR_value", "timepoint", "LR_coord")
  
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
          axis.text.x = element_text(angle=45)) +
    ylab("ALR coord")
  if(is.null(save_name)) {
    show(p)
  } else {
    ggsave(paste0("plots/basset/",level,"/",save_name,".png"), scale=2, width=12, height=2, units="in", dpi=100)
  }
}

# read in and filter full data set at this phylogenetic level
glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)

# cut this down to the desired individual
indiv_data <- subset_samples(data, sname==baboon)

if(vizualization) {
  # visualize the data at this level
  plot_timecourse_phyloseq(indiv_data, save_filename=paste0(baboon,"_phylum_timecourse"), gapped=FALSE, legend=TRUE, legend_level="phylum")
  plot_timecourse_phyloseq(indiv_data, save_filename=paste0(baboon,"_phylum_timecourse_gapped"), gapped=TRUE, legend=TRUE, legend_level="phylum")

  # visualize autocorrelation
  metadata <- read_metadata(data)
  lags <- calc_autocorrelation(data,
                             metadata,
                             lag.max=36,
                             date_diff_units="months",
                             resample=FALSE,
                             use_alr=TRUE,
                             alr_ref=NULL)
  plot_mean_autocorrelation(lags,
                             filename=paste("plots/autocorrelation_36months_GPdiagnostic",sep=""),
                             width=10,
                             height=4)
}

# get observations as differences from baseline in units of days
indiv_metadata <- read_metadata(indiv_data)
baseline_date <- indiv_metadata$collection_date[1]
observations <- sapply(indiv_metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1
Y <- otu_table(indiv_data)@.Data

# chop down to span of interest (if applicable)
if(!is.null(date_lower_limit) & !is.null(date_upper_limit)) {
  # require both be present for now
  min_idx <- min(which(names(observations) >= date_lower_limit))
  max_idx <- max(which(names(observations) <= date_upper_limit))
  # subset dates
  Y_pre <- Y
  Y <- t(Y_pre[min_idx:max_idx,])
  rm(Y_pre)
  # set first observation to t=1
  observations_pre <- observations
  observations <- matrix(observations_pre[min_idx:max_idx], nrow=1) - observations_pre[min_idx] + 1
  rm(observations_pre)
} else {
  # just clean up
  Y <- t(Y)
  dim(observations) <- c(1, length(observations))
}
colnames(Y) <- NULL
rownames(Y) <- NULL

D <- nrow(Y)
N <- ncol(Y)

# abundances <- otu_table(data)@.Data
# alr_abundances <- driver::alr(abundances+0.65, d=alr_ref)
# max_min_diff <- abs(apply(alr_abundances, 2, max) - apply(alr_abundances, 2, min))
# scale sigma by 1/3 since Sigma (posterior taxonomic covariance) contributes to log ratio
# variance *twice*; the idea here is the diagonal of gamma will account for about 1/3 of the
# observed variance in log ratios; this probably isn't quite right but let's see if it's close???
# sigma <- sqrt(median(max_min_diff)) / 3
# scale down by larger size of Gamma relative to Sigma
# sigma <- D/N
sigma <- 1
cat("Using sigma:",sigma,"\n")

# we can make informed(ish) choices for kernel parameters
c <- 0.1 # desired correlation
dd <- 120 # distance in days where we should hit that correlation
rho_se <- sqrt(-dd^2/(2*log(c)))
cat("Using bandwidth (squared exponential):",rho_se,"\n")
# simulate
# x <- 1:730 # distance
# plot(x, (sigma^2)*exp(-(x^2)/(2*rho^2)), type="l", ylim=c(0, 1))

period <- 365
c <- 0.1 # desired correlation
dd <- 180 # distance in days where we should hit that correlation
rho_per <- sqrt(-2*(sin(pi*dd/period)^2)/(log(c)))
cat("Using bandwidth (periodic):",rho_per,"\n")
# simulate
# x <- 1:730 # distance
# lines(x, 0.15*sigma^2*exp(-2*(sin(pi*x/period)^2)/(rho^2)), col="blue")

Gamma <- function(X) 0.15*PER(X, sigma=sigma, rho=rho_per, period=period, jitter=0) + 0.85*SE(X, sigma=sigma, rho=rho_se, jitter=0) + (1e-8)*diag(ncol(X)) # pretty arbitrary

fit_obj <- fit_to_baboon(baboon, Y, observations, Gamma, alr_ref=alr_ref)
predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=100) # interpolates

LR_coords <- c(1,2,3)
plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[1], save_name=paste0(baboon,"_",LR_coords[1]))
plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[2], save_name=paste0(baboon,"_",LR_coords[2]))
plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[3], save_name=paste0(baboon,"_",LR_coords[3]))

Sigma <- fit_obj$fit$Sigma
save(Sigma, file=paste0("subsetted_indiv_data/",level,"/",baboon,"_bassetfit.RData"))




