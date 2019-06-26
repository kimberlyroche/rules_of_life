#library(stray)
#devtools::load_all("/data/mukherjeelab/labraduck")

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

fit_to_baboon <- function(baboon, indiv_data, Gamma, date_lower_limit=NULL, date_upper_limit=NULL, alr_ref=NULL) {
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
    ggsave(paste0("plots/basset/",save_name,".png"), scale=2, width=12, height=3, units="in", dpi=100)
  }
}

# for reference, individuals passable as arguments are:
# "DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI"

if(FALSE) {
args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Usage: Rscript basset_walkthrough.r ACA", call.=FALSE)
}

baboon <- args[1]
}

baboon <- "ACA"

load(paste0("subsetted_indiv_data/",baboon,"_data.RData"))

# parameter settings for GP kernels
# sigma should be approx. the standard deviation between (ALR) samples

if(FALSE) {
glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
alr_ref <- 9 # low-count cohort
# get mean abundances
abundances <- otu_table(data)@.Data
abundance_means <- colMeans(abundances)
alr_mean <- driver::alr(abundance_means, d=alr_ref)
# estimate the variance
total_samples <- nrow(abundances)
alr_abundances <- driver::alr(abundances+0.65, d=alr_ref)
centered_abundances <- t(alr_abundances) - matrix(t(alr_mean), nrow=length(alr_mean), ncol=total_samples)
alr_var <- mean(apply(centered_abundances, 2, function(x) t(x)%*%x))
sigma <- sqrt(alr_var)
}

# we can make informed(ish) choices for kernel parameters
sigma <- 1
c <- 0.1 # desired correlation
dd <- 60 # distance in days where we should hit that correlation
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

alr_ref <- 9

Gamma <- function(X) 0.15*PER(X, sigma=sigma, rho=rho_per, period=period, jitter=0) + 0.85*SE(X, sigma=sigma, rho=rho_se, jitter=0) + (1e-8)*diag(ncol(X)) # pretty arbitrary

fit_obj <- fit_to_baboon(baboon, indiv_data, Gamma, date_lower_limit="2001-10-01", date_upper_limit="2003-11-30", alr_ref=alr_ref)
#fit_obj <- fit_to_baboon(baboon, indiv_data, Gamma)
predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=100) # interpolates

plot_predictions(fit_obj, predict_obj, LR_coord=3, save_name=NULL)

if(FALSE) {
LR_coords <- NULL
# chosen because they give (1) a reference (2) an apparent positive covary-er (3) zero covary-er
if(baboon == "ACA") {
  LR_coords <- c(19, 21, 20)
} else if(baboon == "DUX") {
  LR_coords <- c(19, 21, 20)
} else if(baboon == "LOG") {
  LR_coords <- c(1, 15, 7)
} else if(baboon == "THR") {
  LR_coords <- c(1, 15, 3)
} else if(baboon == "VAI") {
  LR_coords <- c(1, 25, 21)
}

if(!is.null(LR_coords)) {
  plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[1], save_name=paste0(baboon,"_ref"))
  plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[2], save_name=paste0(baboon,"_pos"))
  plot_predictions(fit_obj, predict_obj, LR_coord=LR_coords[3], save_name=paste0(baboon,"_neu"))
}

save(fit_obj, file=paste0("subsetted_indiv_data/",baboon,"_bassetfit.RData"))
}




