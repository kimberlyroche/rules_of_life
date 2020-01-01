library(stray)

source("include/R/general.R")
source("include/R/GP.R")

# periodic kernel
#   X is Q x N as in other kernels
#   rho is bandwidth
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
}

level="family"
date_lower_limit=NULL
date_upper_limit=NULL
alr_ref=38
verbose=TRUE
mean_only=FALSE
max_iter=10
eps_f=NULL
eps_g=NULL

# read in and filter full data set at this phylogenetic level
glom_data <- load_glommed_data(level=level, replicates=TRUE)
subsetted_data <- filter_data(glom_data, level=level, verbose=FALSE)
data <- subset_samples(subsetted_data, sname %in% best_sampled)
metadata <- read_metadata(data)

do_joint <- TRUE
do_independent <- FALSE

if(do_joint) {
  #se_weight=sqrt(2.577)
  #per_weight=sqrt(0.286)
  #wn_weight=0
  se_weight=2
  dd_se=90

  # get observation labels -- space these out by individual
  observations <- metadata$collection_date
  baboons <- unique(metadata$sname)
  individual_start_end <- list()
  for(b in 1:length(baboons)) {
    baboon <- baboons[b]
    idx <- metadata$sname == baboon
    baseline_date <- min(metadata$collection_date[idx])
    temp <- sapply(metadata$collection_date[idx], function(x) round(difftime(x, baseline_date, units="days"))) + 1
    adj_temp <- temp + (b-1)*10000
    observations[idx] <- adj_temp
    individual_start_end[[baboon]] <- c(adj_temp[1], adj_temp[length(adj_temp)])
  }
  Y <- otu_table(data)@.Data
  Y <- t(Y) # D x N samples
  observations <- as.numeric(observations)
  dim(observations) <- c(1, length(observations))

  colnames(Y) <- NULL
  rownames(Y) <- NULL
  D <- nrow(Y)
  N <- ncol(Y)

  # stray uses the D^th element as the ALR reference by default
  # do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }

  if(!file.exists("censor.rds")) {
    censor <- as.logical(rbinom(N, 1, 0.1))
    saveRDS(censor, "censor.rds")
  } else {
    censor <- readRDS("censor.rds")
  }

  Y.train <- Y[,!censor]
  Y.test <- Y[,censor]
  assignments.train <- metadata$sname[!censor]
  assignments.test <- metadata$sname[censor]
  N.train <- ncol(Y.train)
  N.test <- ncol(Y.test)
  observations.train <- observations[,!censor,drop=F]
  observations.test <- observations[,censor,drop=F]
  cat("Censoring",sum(censor),"observations\n")

  # square exponential kernel parameters
  dc <- 0.1 # desired minimum correlation
  se_sigma <- 1
  rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay

  # periodic kernel parameters
  #period <- 365
  #per_sigma <- 1
  #rho_per <- 1

  Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
    #per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
    #wn_weight*WHITENOISE(X, sigma=1, jitter=0) +
    (1e-9)*diag(ncol(X)) # pretty arbitrary

  prior_obj <- default_ALR_prior(D)

  alr_ys <- driver::alr((t(Y.train)+pc))
  # means are per-individual now
  alr_means <- list()
  #alr_means <- matrix(NA, nrow(Y.train)-1, ncol(Y.train))
  individual_lookup <- as.character(max(observations)) # this is essentially the inverse of individual_start_end
                                                       # date index to baboon sname lookup
  for(b in 1:length(baboons)) {
    baboon <- baboons[b]
    idx <- assignments.train == baboon
    temp <- colMeans(alr_ys[idx,])
    alr_means[[baboon]] <- temp
    #alr_means[,idx] <- temp
    individual_lookup[individual_start_end[[baboon]][[1]]:individual_start_end[[baboon]][[2]]] <- baboon
  }

  Theta <- function(X) {
    matrix(0, D-1, ncol(X))
  }

  #Theta <- function(X) {
  #  temp <- matrix(NA, D-1, ncol(X))
  #  for(j in 1:ncol(X)) {
  #    temp[,j] <- alr_means[[individual_lookup[X[1,j]]]]
  #  }
  #  return(temp)
  #}

  if(!file.exists("fit.rds")) {
    ll <- NaN
    while(is.nan(ll)) {
      cat("Fitting model...\n")
      fit <- stray::basset(Y.train, observations.train, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=100)
      ll <- fit$logMarginalLikelihood
    }
    cat("Saving model...\n")
    saveRDS(fit, "fit.rds")
  } else {
    fit <- readRDS("fit.rds")
  }

  # might want to separate these
  X_predict <- c()
  for(baboon in names(individual_start_end)) {
    X_predict <- c(X_predict, (individual_start_end[[baboon]][1]):(individual_start_end[[baboon]][[2]]))
  }
  #X_predict <- X_predict[1:100]
  dim(X_predict) <- c(1, length(X_predict))
  Y.predicted <- predict(fit, X_predict, response="Y", iter=1000)

  baboon.RMSE <- c()
  #idx <- which(X_predict %in% observations.test)
  idx <- sapply(observations.test, function(x) which(X_predict == x))
  for(k in 1:length(idx)) {
    true_counts <- Y.test[,k]
    names(true_counts) <- NULL
    for(m in 1:dim(Y.predicted)[3]) {
      temp <- sqrt(mean((true_counts - Y.predicted[,k,m])^2))
      baboon.RMSE <- c(baboon.RMSE, temp)
    }
  }
  all.RMSE <- rbind(all.RMSE, data.frame(RMSE=baboon.RMSE, baboon=rep(baboon, length(baboon.RMSE))))

  p <- ggplot(all.RMSE, aes(x=RMSE, color=baboon)) +
      geom_density()
    ggsave("test_joint.png", plot=p, dpi=200, units="in", width=15, height=10)

}

if(do_independent) {

  #se_weight=sqrt(2.577)
  #per_weight=sqrt(0.286)
  #wn_weight=0
  se_weight=sqrt(2.863)
  dd_se=90

  all.RMSE <- data.frame(RMSE=c(), baboon=c())
  for(baboon in best_sampled) {
    cat("Baboon:",baboon,"\n")
    # global var hack still necessary for phyloseq subset_samples (I think)
    baboon <<- baboon
    # cut this down to the desired individual
    indiv_data <- subset_samples(data, sname==baboon)

    # get observations as differences from baseline in units of days
    indiv_metadata <- read_metadata(indiv_data)
    baseline_date <- indiv_metadata$collection_date[1]
    observations <- sapply(indiv_metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1
    Y <- otu_table(indiv_data)@.Data
    Y <- t(Y)
    dim(observations) <- c(1, length(observations))

    colnames(Y) <- NULL
    rownames(Y) <- NULL
    D <- nrow(Y)
    N <- ncol(Y)

    # stray uses the D^th element as the ALR reference by default
    # do some row shuffling in Y to put the reference at the end
    if(!is.null(alr_ref)) {
      Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
    }

    censor <- as.logical(rbinom(N, 1, 0.1))

    Y.train <- Y[,!censor]
    Y.test <- Y[,censor]
    N.train <- ncol(Y.train)
    N.test <- ncol(Y.test)
    observations.train <- observations[,!censor,drop=F]
    observations.test <- observations[,censor,drop=F]
    cat("Censoring",sum(censor),"observations\n")
      
    # square exponential kernel parameters
    dc <- 0.1 # desired minimum correlation
    se_sigma <- 1
    rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay

    # periodic kernel parameters
    #period <- 365
    #per_sigma <- 1
    #rho_per <- 1

    Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
      #per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
      #wn_weight*WHITENOISE(X, sigma=1, jitter=0) +
      (1e-8)*diag(ncol(X)) # pretty arbitrary

    prior_obj <- default_ALR_prior(D)

    alr_ys <- driver::alr((t(Y.train)+pc))
    alr_means <- colMeans(alr_ys)
    Theta <- function(X) matrix(alr_means, D-1, ncol(X))

    ll <- NaN
    while(is.nan(ll)) {
      cat("Fitting model...\n")
      fit <- stray::basset(Y.train, observations.train, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=1000)
      ll <- fit$logMarginalLikelihood
    }

    # note: predictions blow up (specifically Gamma_ooIou blows up) if I don't predict with some
    # known observations; not sure if this is a bug or a conceptual problem with what I'm trying
    # to do but for now, we'll predict a whole trajectory (observed and censored samples) and
    # pick out the previously censored

    X_predict <- t(1:max(observations))
    Y.predicted <- predict(fit, X_predict, response="Y", iter=1000)

    baboon.RMSE <- c()
    #idx <- which(X_predict %in% observations.test)
    idx <- sapply(observations.test, function(x) which(X_predict == x))
    for(k in 1:length(idx)) {
      true_counts <- Y.test[,k]
      names(true_counts) <- NULL
      for(m in 1:dim(Y.predicted)[3]) {
        temp <- sqrt(mean((true_counts - Y.predicted[,k,m])^2))
        baboon.RMSE <- c(baboon.RMSE, temp)
      }
    }
    all.RMSE <- rbind(all.RMSE, data.frame(RMSE=baboon.RMSE, baboon=rep(baboon, length(baboon.RMSE))))
  }

  p <- ggplot(all.RMSE, aes(x=RMSE, color=baboon)) +
      geom_density()
    ggsave("test.png", plot=p, dpi=200, units="in", width=15, height=10)

}

