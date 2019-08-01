library(stray)
#devtools::load_all("/data/mukherjeelab/Mongrel/stray")
#devtools::load_all("/data/mukherjeelab/Mongrel/labraduck")
library(driver)

# for reference, individuals passable as arguments are:
# "DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI"

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  # arguments are+
  #   baboon sname
  #   agglomeration level
  #   SE scale
  #   PER scale
  #   WN scale
  #   days decay for SE kernel
  #   ALR reference taxon index
  #   plot save name append string
  stop("Usage: Rscript 01_fit_basset.R ECH family 1.98 0.22 0 90", call.=FALSE)
}
baboon <- args[1]
level <- args[2]
if(length(args) >= 3) {
  se_weight <- as.numeric(args[3])
  per_weight <- se_weight*as.numeric(args[4])
  wn_weight <- se_weight*as.numeric(args[5])
} else {
  se_weight <- 2
  per_weight <- 0.25
  wn_weight <- 0
}
if(length(args) == 6) {
  dd_se <- as.numeric(args[6])
} else {
  dd_se <- 90
}
if(length(args) == 7) {
  alr_ref <- as.numeric(args[7])
} else {
  alr_ref <- NULL
}
if(length(args) >= 8) {
  save_append <- paste0("_",args[8])
} else {
  save_append <- ""
}

cat(paste0("Running stray::basset with with parameters:\n\tbaboon=",baboon,"\n\tlevel=",level,"\n\tSE kernel weight=",se_weight,"\n\tPER kernel weight=",per_weight,"\n\tdd_se=",dd_se,"\n\talr_ref=",alr_ref,"\n"))

# testing
#date_lower_limit <- "2001-10-01"
#date_upper_limit <- "2003-01-01"
#date_lower_limit <- "2004-05-24" # this includes a big gap if run on ACA
#date_upper_limit <- "2008-10-25" # useful to check effect of prior mean
date_lower_limit <- NULL
date_upper_limit <- NULL

# generate plots for these (basically random) individuals, for diagnostics
plot_these <- c("POW", "DUI", "COO", "YAI", "ACA", "ZIZ")

library(phyloseq)
library(dplyr)
library(driver)
library(tidyverse)

source("include.R")

get_predictions <- function(X, fit, n_samples=100){
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
  Xi <- GG%*%(diag(D)*1)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
  
  alr_ys <- driver::alr((t(Y)+0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
  fit <- stray::basset(Y, observations, upsilon, Theta, Gamma, Xi)
  return(list(Y=Y, alr_ys=alr_ys, X=observations, fit=fit))
}

plot_predictions <- function(fit_obj, predict_obj, LR_coord=1, save_name=NULL) {
  observations <- fit_obj$X
  lr_tidy <- gather_array(fit_obj$alr_ys, "LR_value", "timepoint", "LR_coord")  
  
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
    ylab("LR coord")
  if(is.null(save_name)) {
    show(p)
  } else {
    ggsave(paste0("plots/basset/",level,"/",save_name,".png"), scale=2, width=12, height=2, units="in", dpi=100)
  }
}

# read in and filter full data set at this phylogenetic level
glom_data <- load_glommed_data(level=level, replicates=TRUE)
# previous thresholding: data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
# this thresholding was chosen on the basis of retaining 98.7% of counts
# it retains taxa (1) with taxonomic identification to at least level=level and (2) which are present
#     at least a 5-count in at least 1/3 of samples
data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.33, collapse_level=level, verbose=TRUE)

# cut this down to the desired individual
indiv_data <- subset_samples(data, sname==baboon)

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

#dd_se <- 30
dc <- 0.1 # desired SE minimum correlation
se_sigma <- 1
rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate

period <- 365
per_sigma <- 1
rho_per <- 1

Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
                     per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
                     wn_weight*WHITENOISE(X, sigma=1, jitter=0) +
                     (1e-8)*diag(ncol(X)) # pretty arbitrary

fit_obj <- fit_to_baboon(baboon, Y, observations, Gamma, alr_ref=alr_ref)
save(fit_obj, file="test.RData")

predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=100)

LR_coords <- c()
if(level == "phylum") {
  LR_coords <- c(1,2,3,7,8)
} else if(level == "family") {
  LR_coords <- c(1,2,3,20,21)
}
for(coord in LR_coords) {
  if(baboon %in% plot_these) {
    plot_predictions(fit_obj, predict_obj, LR_coord=coord, save_name=paste0(baboon,"_",coord,save_append))
  }
}

subset_sz <- 100

if(FALSE) {
Eta <- fit_obj$fit$Eta[,,1:subset_sz]
Lambda <- fit_obj$fit$Lambda[,,1:subset_sz]
Sigma <- fit_obj$fit$Sigma[,,1:subset_sz]
save(Eta, Lambda, Sigma, file="test_ALR.RData")
}

V <- driver::create_default_ilr_base(ncategories(fit_obj$fit))
fit.ilr <- to_ilr(fit_obj$fit, V)
Eta <- fit.ilr$Eta[,,1:subset_sz]
Lambda <- fit.ilr$Lambda[,,1:subset_sz]
Sigma <- fit.ilr$Sigma[,,1:subset_sz]
#save(Eta, Lambda, Sigma, file="test_ILR.RData")

if(FALSE) {
fit.clr <- to_clr(fit_obj$fit)
Eta <- fit.clr$Eta[,,1:subset_sz]
Lambda <- fit.clr$Lambda[,,1:subset_sz]
Sigma <- fit.clr$Sigma[,,1:subset_sz]
save(Eta, Lambda, Sigma, file="test_CLR.RData")
}

save(V, Eta, Lambda, Sigma, file=paste0("subsetted_indiv_data/",level,"/",baboon,"_bassetfit",save_append,".RData"))
