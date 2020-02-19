# this file fits the model over all individuals together to estimate a single covariance
# matrix associated with all taxa across all hosts

relative_path <- ".."

# usage: Rscript joint.R {optional: k, where k specifies the use of the first k best sampled hosts}

args <- commandArgs(trailingOnly=TRUE)
individual_limit <- 10
if(length(args) >= 1) {
  individual_limit <- as.numeric(args[1])
}
cat("Using",individual_limit,"individuals...\n")

source(file.path(relative_path,"include/R/general.R"))
source(file.path(relative_path,"include/R/GP.R"))

do_save <- TRUE
predict_level <- 1 # 1 = only predict on hold-outs; calculate RMSE distribution
                   # 2 = predict all observations (train and test); render prediction error plots
                   # 3 = predict all observations with day-level interpolation; unused for now

level <- "family"
alr_ref <- 38

glom_data <- load_glommed_data(level=level, replicates=TRUE)
filtered_data <- filter_data(glom_data, level=level, verbose=FALSE)
if(individual_limit >= 10) {
  indiv_list <- best_sampled
} else {
  indiv_list <- best_sampled[1:individual_limit]
}

data <- subset_samples(filtered_data, sname %in% indiv_list)
metadata <- read_metadata(data)

N <- length(metadata$collection_date)

# get time indices for each individual (space these out by margin)
# record start and end indices for each individual
observations <- matrix(NA, 2, N)
observations[1,] <- metadata$collection_date
baboons <- unique(metadata$sname)
individual_start_end <- list()
margin <- 10000
for(b in 1:length(baboons)) {
  # for each baboon
  baboon <- baboons[b]
  idx <- metadata$sname == baboon
  baseline_date <- min(metadata$collection_date[idx])
  # get time indices associated with collection dates (1-indexed)
  temp <- sapply(metadata$collection_date[idx], function(x) round(difftime(x, baseline_date, units="days"))) + 1
  # add margin
  adj_temp <- temp + (b-1)*margin
  observations[1,idx] <- adj_temp
  observations[2,idx] <- b
  # record the start and end indices for this individual
  individual_start_end[[b]] <- c(adj_temp[1], adj_temp[length(adj_temp)])
}
Y <- otu_table(data)@.Data
Y <- t(Y) # dimensions are: D taxa x N samples
observations <- apply(observations, c(1,2), as.numeric)

if(do_save) {
  saveRDS(individual_start_end, file.path(relative_path,output_dir,"joint_boundaries.rds"))
}

colnames(Y) <- NULL
rownames(Y) <- NULL
D <- nrow(Y)

# stray uses the D^th element as the ALR reference by default
# do some row shuffling in Y to put the reference at the end
if(!is.null(alr_ref)) {
  Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
}

censor <- NULL
filename <- file.path(relative_path,output_dir,"joint_censor.rds")
if(do_save & file.exists(filename)) {
  censor <- readRDS(filename)
}
if(is.null(censor)) {
  censor <- as.logical(rbinom(N, 1, 0.1))
  if(do_save) {
    saveRDS(censor, filename)
  }  
}

# get train/test split
Y.train <- Y[,!censor]
Y.test <- Y[,censor]
N.train <- ncol(Y.train)
N.test <- ncol(Y.test)
observations.train <- observations[,!censor,drop=F]
observations.test <- observations[,censor,drop=F]
cat("Censoring",sum(censor),"observations\n")

# set up kernel stuff

dd_se <- 90
dc <- 0.1 # desired minimum correlation
rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
se_sigma <- 2

Gamma <- function(X) {
  SE(X[1,,drop=F], sigma=se_sigma, rho=1, jitter=1e-10)
}

prior_obj <- default_ALR_prior(D)

# calculate the mean ALR vector for each individual
alr_ys <- t(driver::alr((t(Y.train)+pc))) # dimensions are: D-1 taxa x N samples
# index individuals by their index in the {baboon} vector
alr_means <- list()
for(b in 1:length(baboons)) {
  baboon <- baboons[b]
  idx <- observations.train[2,] == b
  temp <- rowMeans(alr_ys[,idx])
  alr_means[[b]] <- temp
}

Theta <- function(X) {
  temp <- matrix(NA, D-1, ncol(X))
  for(j in 1:ncol(X)) {
    temp[,j] <- alr_means[[X[2,j]]]
  }
  return(temp)
}

fit <- NULL
filename <- file.path(relative_path,output_dir,"joint_fit.rds")
if(do_save & file.exists(filename)) {
  fit <- readRDS(filename)
}
if(is.null(fit)) {
  cat("Fitting model...\n")
  fit <- stray::basset(Y.train, observations.train, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=100)
  ll <- fit$logMarginalLikelihood
  cat("\tLog likelihood:",ll,"\n")
  if(do_save) {
    saveRDS(fit, filename)
  }  
}

Y.predicted <- NULL
filename <- file.path(relative_path,output_dir,"joint_predict.rds")
#if(do_save & file.exists(filename)) {
#  Y.predicted <- readRDS(filename)
#}
#if(is.null(Y.predicted)) {
  cat("Predicting from model...\n")
  if(predict_level == 3) {
    N.predict <- sum(sapply(individual_start_end, function(x) x[2]-x[1]+1))
    X_predict <- matrix(NA, 2, N.predict)
    offset <- 1
    for(b in 1:length(baboons)) {
      b_range <- individual_start_end[[b]][1]:individual_start_end[[b]][2]
      b_range_len <- length(b_range)
      X_predict[1,offset:(offset+b_range_len-1)] <- b_range
      X_predict[2,offset:(offset+b_range_len-1)] <- b
      offset <- offset + b_range_len
    }
    Y.predicted <- predict(fit, X_predict, response="Y")
  } else if(predict_level == 2) {
    Y.predicted <- predict(fit, observations, response="Y")
  } else if(predict_level == 1) {
    Y.predicted <- predict(fit, observations.test, response="Y")
  }
  if(do_save) {
    saveRDS(Y.predicted, filename)
  }
#}

if(predict_level == 1) {
  # indexing here only works for predictions on holdouts right now
  # for full predictions, we'll need to get the indices in X_predict that correspond to
  #   the holdouts we want to test
  baboon.RMSE <- data.frame(RMSE=c())
  idx <- observations.test[1,]
  for(k in 1:length(idx)) {
    true_counts <- Y.test[,k]
    names(true_counts) <- NULL
    for(m in 1:dim(Y.predicted)[3]) {
      # for each sample prediction
      temp <- sqrt(mean((true_counts - Y.predicted[,k,m])^2))
      baboon.RMSE <- rbind(baboon.RMSE, data.frame(RMSE=temp))
    }
  }
  p <- ggplot(baboon.RMSE, aes(x=RMSE)) +
      geom_density()
  ggsave(file.path(relative_path,output_dir,"evaljoint_errorplot_joint.png"), plot=p, dpi=200, units="in", width=15, height=10)
}

if(predict_level == 2) {
  # collect the ACTUAL counts
  ground_truth <- data.frame(observation=c(), taxon=c(), individual=c(), count=c(), holdout=c())
  for(d in 1:D) {
    ground_truth <- rbind(ground_truth, data.frame(observation=observations[1,],
                                                   taxon=d,
                                                   individual=observations[2,],
                                                   count=Y[d,],
                                                   holdout=censor))
  }
  ground_truth$sample <- 0
  ground_truth$count <- ground_truth$count + 0.5
  ground_truth <- ground_truth %>%
    mutate(log_count = log10(count))

  # collect the PREDICTED counts
  predicted_tidy <- gather_array(Y.predicted, "count", "taxon", "timepoint", "sample")
  # convert linear time steps to observations indices
  map <- data.frame(timepoint=1:ncol(observations), observation=observations[1,])
  predicted_tidy <- merge(predicted_tidy, map, by="timepoint")
  predicted_tidy <- predicted_tidy[,!(names(predicted_tidy) %in% c("timepoint"))]
  # add holdout flag
  map <- data.frame(observation=observations[1,], holdout=censor)
  predicted_tidy <- merge(predicted_tidy, map, by="observation")
  # add individual label
  map <- data.frame(observation=observations[1,], individual=observations[2,])
  predicted_tidy <- merge(predicted_tidy, map, by="observation")
  # convert counts to log(counts + pseudocount)
  predicted_tidy$count <- predicted_tidy$count + 0.5
  predicted_tidy <- predicted_tidy %>%
    mutate(log_count=log10(count))

  predicted_quantiles <- predicted_tidy %>%
    group_by(individual, taxon, observation) %>%
    summarise(p25 = quantile(log_count, prob=0.25),
              mean = mean(log_count),
              p75 = quantile(log_count, prob=0.75)) %>%
  ungroup()
  predicted_quantiles <- as.data.frame(predicted_quantiles)

  p <- ggplot() +
    geom_errorbar(data=predicted_quantiles[predicted_quantiles$taxon == 1,],
              aes(x=observation, ymin=p25, ymax=p75), width=0.5, color="#666666") +
    geom_point(data=ground_truth[ground_truth$taxon == 1,],
              aes(x=observation, y=log_count, color=holdout), size=1) +
    facet_wrap(~individual, ncol=2, scales="free_x") +
    theme(legend.position = "none")
  ggsave(file.path(relative_path,output_dir,"evaljoint_RMSE.png"), p, units="in", dpi=100, height=10, width=20)
}

