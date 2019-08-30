source("include/R/general.R")
source("include/R/metagenomics.R")
sourceCpp("include/cpp/Riemann_dist.cpp")

GP_plot_dir <- paste0(plot_dir,"basset/")

# ====================================================================================================================
# ADDITIONAL GP KERNELS
# ====================================================================================================================

# periodic kernel
#   X is Q x N as in other kernels
#   rho is bandwidth
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
}

# whitenoise kernel
WHITENOISE <- function(X, sigma=1, jitter=1e-10) {
  dist <- as.matrix(dist(t(X)))
  G <- diag(ncol(dist))*sigma^2 + jitter*diag(ncol(dist))
  return(G)
}

# ====================================================================================================================
# STRAY WRAPPERS
# ====================================================================================================================

# fit Gaussian process to a single baboon series using stray::basset
fit_GP <- function(baboon, level, se_weight, per_weight, wn_weight, dd_se, save_append="",
                   date_lower_limit=NULL, date_upper_limit=NULL, alr_ref=NULL, verbose=TRUE, mean_only=FALSE) {
  if(verbose) {
    cat(paste0("Fitting stray::basset with with parameters:\n",
             "\tbaboon=",baboon,"\n",
             "\tlevel=",level,"\n",
             "\tSE kernel weight=",se_weight,"\n",
             "\tPER kernel weight=",per_weight,"\n",
             "\tWN kernel weight=",wn_weight,"\n",
             "\tdd_se=",dd_se,"\n"))
  }

  # read in and filter full data set at this phylogenetic level
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, collapse_level=level, verbose=FALSE)
  
  # global var hack still necessary for phyloseq subset_samples (I think)
  baboon <<- baboon
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

  alr_ref <- pick_alr_ref(D)
  
  # square exponential kernel parameters
  dc <- 0.1 # desired minimum correlation
  se_sigma <- 1
  rho_se <- sqrt(-dd_se^2/(2*log(dc))) # back calculate the decay
  
  # periodic kernel parameters
  period <- 365
  per_sigma <- 1
  rho_per <- 1
  
  Gamma <- function(X) se_weight*SE(X, sigma=se_sigma, rho=rho_se, jitter=0) +
    per_weight*PER(X, sigma=per_sigma, rho=rho_per, period=period, jitter=0) +
    wn_weight*WHITENOISE(X, sigma=1, jitter=0) +
    (1e-10)*diag(ncol(X)) # pretty arbitrary
  
  D <- nrow(Y)
  N <- ncol(Y)
  
  # stray uses the D^th element as the ALR reference by default
  # do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }
  
  prior_obj <- default_ALR_prior(D)
  
  alr_ys <- driver::alr((t(Y)+pc))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
  if(mean_only) {
    fit <- stray::basset(Y, observations, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi, n_samples=0, ret_mean=TRUE)
  } else {
    fit <- stray::basset(Y, observations, prior_obj$upsilon, Theta, Gamma, prior_obj$Xi)
  }
  fit_obj <- list(Y=Y, alr_ys=alr_ys, X=observations, fit=fit)
  
  # dumb as hell but for later prediction, these need to be loaded into the workspace
  # for use by Gamma; fix this eventually
  fit_obj$kernelparams$se_weight <- se_weight
  fit_obj$kernelparams$se_sigma <- se_sigma
  fit_obj$kernelparams$rho_se <- rho_se
  fit_obj$kernelparams$per_weight <- per_weight
  fit_obj$kernelparams$per_sigma <- per_sigma
  fit_obj$kernelparams$rho_per <- rho_per
  fit_obj$kernelparams$period <- period
  fit_obj$kernelparams$wn_weight <- wn_weight
  
  # chop down for size savings since /data/mukherjeelab is full as shit
  # can't seem to pass desired sample number to stray with any effect; debug this eventually
  fit_obj$fit$iter <- 100
  if(mean_only) {
    fit_obj$fit$iter <- 1
  }
  fit_obj$fit$Eta <- fit_obj$fit$Eta[,,1:fit_obj$fit$iter]
  fit_obj$fit$Lambda <- fit_obj$fit$Lambda[,,1:fit_obj$fit$iter]
  fit_obj$fit$Sigma <- fit_obj$fit$Sigma[,,1:fit_obj$fit$iter]

  if(mean_only) {
    # just return this small, summary of the GP fit
    return(fit_obj)
  }

  # otherwise, save the full posterior sample set
  saveRDS(fit_obj, paste0(model_dir,level,"/",baboon,"_bassetfit",save_append,".rds"))
}

# predict from fitted stray::basset model; predictions are in CLR by default
#   X is the observation vector in days (e.g. 1,2,17,21,22...); this will be interpolated
#   fit is a strayfit/bassetfit object
get_predictions <- function(X, fit, n_samples=100) {
  cat("Predicting from 1 to",max(X),"\n")
  X_predict <- t(1:(max(X)))
  fit.clr <- to_clr(fit)
  predicted <- predict(fit.clr, X_predict, iter=n_samples) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

# ====================================================================================================================
# FITTED MODEL LIST ACCESSORS
# ====================================================================================================================

# get filename list for all fitted basset models
get_fitted_modellist <- function(level="family") {
  pattern_str <- "*_bassetfit.rds"
  regexpr_str <- "_bassetfit.rds"
  fitted_models <- list.files(path=paste0(model_dir,level), pattern=pattern_str, full.names=TRUE, recursive=FALSE)
  return(list(fitted_models=fitted_models,
              pattern_str=pattern_str,
              regexpr_str=regexpr_str))
}

# get filename list for all fitted basset models and do some extra result parsing
get_fitted_modellist_details <- function(level="family") {
  fitted_models <- get_fitted_modellist(level)
  individuals <- sapply(fitted_models$fitted_models, function(x) { idx <- regexpr(fitted_models$regexpr_str, x); return(substr(x, idx-3, idx-1)) } )
  names(individuals) <- NULL
  # parse dimensions
  fit_obj <- readRDS(paste0(model_dir,level,"/",individuals[1],fitted_models$regexpr_str))
  return(list(individuals=individuals,
              pattern_str=fitted_models$pattern_str,
              regexpr_str=fitted_models$regexpr_str,
              model_list=fitted_models$fitted_models,
              N=fit_obj$fit$N,
              D=fit_obj$fit$D,
              n_samples=fit_obj$fit$iter))
}

# ====================================================================================================================
# POSTERIOR EMBEDDING
# ====================================================================================================================

calc_posterior_distances <- function(level, which_measure="Sigma", mean_distance=FALSE) {
  # grab all fitted models
  indiv_obj <- get_fitted_modellist_details(level=level)
  individuals <- indiv_obj$individuals
  pattern_str <- indiv_obj$pattern_str
  regexpr_str <- indiv_obj$regexpr_str
  P <- indiv_obj$D - 1 # ILR, else P = D-1
  N <- indiv_obj$N
  n_samples <- indiv_obj$n_samples
  n_indiv <- length(individuals)
  if(which_measure == "Sigma") {
    all_samples <- matrix(NA, P, P*n_samples*n_indiv)
  } else {
    # we'll use per-sample average to manage individuals having different N so each individual's
    # posterior will be summarized as one mean vector
    all_samples <- matrix(NA, n_indiv*n_samples, P)
  }
  indiv_labels <- c()
  for(i in 1:n_indiv) {
    fit <- readRDS(indiv_obj$model_list[i])$fit
    # to ILR
    V <- driver::create_default_ilr_base(ncategories(fit))
    fit.ilr <- to_ilr(fit, V)
    Lambda <- fit.ilr$Lambda
    Sigma <- fit.ilr$Sigma
    # to CLR
    # fit.clr <- to_clr(fit)
    # Lambda <- fit.clr$Lambda
    # Sigma <- fit.clr$Sigma
    if(which_measure == "Sigma") {
      Sigma <- Sigma[,,1:n_samples]
      all_samples[,((i-1)*P*n_samples+1):(i*P*n_samples)] <- Sigma
      indiv_labels <- c(indiv_labels, rep(individuals[i], n_samples))
    } else {
      collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) })) # 100 x P
      all_samples[((i-1)*n_samples+1):(i*n_samples),] <- collLambda
      indiv_labels <- c(indiv_labels, rep(individuals[i], n_samples))
    }
  }
  
  if(which_measure == "Sigma") {
    distance_mat <- Riemann_dist_samples(all_samples, n_indiv, n_samples)
  } else {
    distance_mat <- dist(all_samples)
  }

  return(distance_mat)
}

# embed posterior samples of (which_measure) using MDS and the appropriate distance metric
#   which_measure: Sigma | Lambda
embed_posteriors <- function(level, which_measure="Sigma") {
  distance_mat <- calc_posterior_distances(level, which_measure=which_measure)

  cat("Embedding posterior samples...\n")  
  fit <- cmdscale(distance_mat, eig=TRUE, k=6) # k is the number of dim
  cat("\tEigenvalue #1:",fit$eig[1],"\n")
  cat("\tEigenvalue #2:",fit$eig[2],"\n")
  cat("\tEigenvalue #3:",fit$eig[3],"\n")
  cat("\tEigenvalue #4:",fit$eig[4],"\n")
  cat("\tEigenvalue #5:",fit$eig[5],"\n")
  cat("\tEigenvalue #6:",fit$eig[6],"\n")
  cat("\tEigenvalue #7:",fit$eig[7],"\n")
  
  # save first 6 coordinates (arbitrarily)
  df <- data.frame(x1=fit$points[,1], x2=fit$points[,2],
                   x3=fit$points[,3], x4=fit$points[,4],
                   x5=fit$points[,5], x6=fit$points[,6],
                   labels=indiv_labels)
  saveRDS(df, paste0(GP_plot_dir,level,"/",which_measure,"_ordination.rds"))
  
  # centroids are useful for labeling plots
  df_centroids <- df %>%
    group_by(labels) %>%
    summarise(mean_x1=mean(x1), mean_x2=mean(x2),
              mean_x3=mean(x3), mean_x4=mean(x4),
              mean_x5=mean(x5), mean_x6=mean(x6))
  saveRDS(df_centroids, paste0(GP_plot_dir,level,"/",which_measure,"_ordination_centroids.rds"))
}

# ====================================================================================================================
# EMBEDDING VISUALIZATION FNS.
# ====================================================================================================================

load_outcomes <- function() {
  outcomes <- read.csv(paste0(data_dir,"fitness/IndividualTraits_ForKim.csv"), header=TRUE)
  # note: may need to change this filtration eventually
  outcomes <- outcomes[outcomes$sname %in% over_50,]
  # filter to NA-less measures
  outcomes <- outcomes[,apply(outcomes, 2, function(x) sum(is.na(x))==0)]
  return(outcomes) # indexed by sname column
}

# get annotation labels
#   df are the coordinates from the ordination
#   data is the full ABRP phyloseq object
#   individuals is a list of snames
#   annotation is the label to grab
get_other_labels <- function(df, data, individuals, annotation="group") {
  labels <- numeric(nrow(df))
  names(labels) <- df$label
  if(annotation == "group") {
    metadata <- sample_data(data)
    primary_group <- metadata %>%
      select(c("sname", "collection_date", "grp")) %>%
      filter(sname %in% individuals) %>% 
      group_by(sname, grp) %>%
      tally() %>%
      slice(which.max(n))
    for(indiv in individuals) {
      labels[df$label == indiv] <- primary_group[primary_group$sname == indiv,]$grp[[1]]
    }
  }
  if(annotation == "matgroup") {
    metadata <- sample_data(data)
    for(indiv in individuals) {
      labels[df$label == indiv] <- metadata[metadata$sname == indiv,]$matgrp[[1]]
    }
  }
  if(annotation == "counts" | annotation == "density") {
    for(indiv in individuals) {
      cat("Parsing individual",indiv,"...\n")
      indiv <<- indiv # needs to be global (bug)
      indiv_subset <- subset_samples(data, sname==indiv)
      sample_count <- phyloseq::nsamples(indiv_subset)
      labels[df$label == indiv] <- round(sample_count, -1) # discretize
    }
  }
  if(annotation == "density") {
    metadata <- sample_data(data)
    md_subset <- metadata %>%
      select(c("sname", "collection_date")) %>%
      filter(sname %in% individuals) %>%
      group_by(sname) %>%
      summarize(delta=difftime(max(collection_date), min(collection_date), units="days"))
    for(indiv in individuals) {
      labels[names(labels) == indiv] <- labels[names(labels) == indiv]/md_subset[md_subset$sname == indiv,]$delta[[1]]
    }
    labels <- round(labels*100)
  }
  if(annotation %in% c("mom", "dad")) {
    outcomes <- load_outcomes()
    for(indiv in individuals) {
      if(annotation == "mom") { label <- as.character(outcomes[outcomes$sname == indiv,]$mom) }
      else { label <- as.character(outcomes[outcomes$sname == indiv,]$dad) }
      if(label == "") { label <- NA }
      labels[names(labels) == indiv] <- label
    }
  }
  if(annotation %in% c("momrank", "drought", "largegroup", "momdied", "competingsib", "earlyadversity")) {
    outcomes <- load_outcomes()
    for(indiv in individuals) {
      if(annotation == "momrank") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$mom_lowQuartRank }
      if(annotation == "drought") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$bornInDrought }
      if(annotation == "largegroup") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$born_largeGroup }
      if(annotation == "momdied") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$mom_died }
      if(annotation == "competingsib") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$has_CompetingSib }
      if(annotation == "earlyadversity") { labels[names(labels) == indiv] <- outcomes[outcomes$sname == indiv,]$EarlyAdversityScore }
    }
    labels <- as.factor(labels)
  }
  if(annotation %in% c("birthrate_all", "birthrate_surviving")) {
    outcomes <- load_outcomes()
    years <- sapply(as.vector(outcomes$birth_date), function(x) {
      year <- as.numeric(strsplit(x, "/")[[1]][3]);
      if(year < 10) { year <- year + 2000 }
      else { year <- year + 1900 }
    })
    names(years) <- outcomes$sname

    metadata <- sample_data(data)

    for(indiv in individuals) {
      if(outcomes[outcomes$sname==indiv,]$sex == "F") {
        years_obs_cont <- difftime(max(metadata[metadata$sname==indiv,]$collection_date), min(metadata[metadata$sname==indiv,]$collection_date), units="weeks")/52
        if(years_obs_cont >= 1) {
          if(annotation == "birthrate_all") {
            births <- outcomes[outcomes$sname==indiv,]$num_live_births_RAW
          } else {
            births <- outcomes[outcomes$sname==indiv,]$num_surv_births_RAW
          }
          score <- births/as.numeric(years_obs_cont)
          labels[names(labels) == indiv] <- round(score,1)
        } else {
          labels[names(labels) == indiv] <- NA
        }
      } else {
        labels[names(labels) == indiv] <- NA
      }
    }
    labels <- as.factor(labels)
  }

  df2 <- data.frame(x1=df$x1, x2=df$x2,
                    x3=df$x3, x4=df$x4,
                    x5=df$x5, x6=df$x6,
                    labels=as.factor(labels))
  return(df2)
}

# plot 2D embedding using desired coordinate pair as (x,y)
#   df are the coordinates from the ordination
#   df_centroids are the individual centroids (optional)
#   axis1 is the x-axis surrogate
#   axis2 is the y-axis surrogate
#   label_type is the annotation/label to grab
plot_axes <- function(df, df_centroids=NULL, axis1="x1", axis2="x2", label_type="individual", legend=TRUE) {
  p <- ggplot() + geom_point(data=df, aes_string(x=axis1, y=axis2, color="labels"))
  if(label_type == "individual") {
    # label the centroids directly
    p <- p + geom_text(data=df_centroids, aes_string(x=paste0("mean_",axis1), y=paste0("mean_",axis2), label="labels"),
                       color="black", fontface="bold")
  }
  if(!legend | label_type == "individual") {
    p <- p + theme(legend.position='none')
  }
  plot_save_name <- paste0(which_measure,"_ordination_",label_type,"_",axis1,axis2,".png")
  img_width <- 4
  if(legend) {
    img_width <- 4.5
  }
  ggsave(paste0(GP_plot_dir,level,"/",plot_save_name), plot=p, scale=2,
         width=img_width, height=4, units="in", dpi=100)
}

# wrapper to plot ordination with 
# allowable (label_types) values are: individual, group, counts, density
plot_ordination <- function(level, which_measure, label_type, legend=TRUE) {
  df <- readRDS(paste0(GP_plot_dir,level,"/",which_measure,"_ordination.rds"))
  df_centroids <- readRDS(paste0(GP_plot_dir,level,"/",which_measure,"_ordination_centroids.rds"))
  if(label_type != "individual") {
    glom_data <- load_glommed_data(level=level, replicates=TRUE)
    df <- get_other_labels(df, glom_data, unique(df$label), annotation=label_type)
  }
  plot_axes(df, df_centroids, "x1", "x2", label_type, legend=legend)
  plot_axes(df, df_centroids, "x2", "x3", label_type, legend=legend)
  plot_axes(df, df_centroids, "x3", "x4", label_type, legend=legend)
  plot_axes(df, df_centroids, "x4", "x5", label_type, legend=legend)
  plot_axes(df, df_centroids, "x5", "x6", label_type, legend=legend)
}

# ====================================================================================================================
# EMBEDDING VISUALIZATION CONT'D. - AXES EXTREMA
# ====================================================================================================================

# plot proportional/bar plots for most extreme individuals in chosen coordinate of BASELINE embedding;
# by default, truncate taxa to those above a minimum proportional abundance so we're not plotting slivers
#   coordinate is 1...K
#   level is taxonomic level (e.g. family)
#   no_indiv are # individuals at each extreme to visualize
plot_extreme_Lambda <- function(coordinate, level, no_indiv=10, save_filename="extreme_Lambda") {
  df_centroids <- readRDS(paste0(GP_plot_dir,level,"/Lambda_ordination_centroids.rds"))
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))
  
  # first pass: get non-tiny proportions to retain
  # if we don't exclude these, the legend will be unparsable
  allLambda <- NULL
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"(1)\n")
    fit_obj <- readRDS(paste0(model_dir,level,"/",baboon,"_bassetfit.rds"))
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
  
  # stupid feature naming is for the benefit of the existing function plot_timecourse_metagenomics
  # need to generalize it
  df <- data.frame(sample=c(), enzyme=c(), proportion=c())
  
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"(2)\n")
    fit_obj <- readRDS(paste0(model_dir,level,"/",baboon,"_bassetfit.rds"))
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
  plot_timecourse_metagenomics(df, save_filename=paste0(GP_plot_dir,level,"/",save_filename))
}

# plot the diagonal of the element-wise mean covariance across individuals in a single heatmap
# the idea here is to compare at a glance the total variation for each log ratio across all individuals
#   level is taxonomic level (e.g. family)
#   lrtransform is clr, alr, ilr
plot_diag_Sigma <- function(level, lrtransform="clr") {
  indiv_obj <- get_fitted_modellist_details(level=level)
  individuals <- indiv_obj$individuals
  pattern_str <- indiv_obj$pattern_str
  regexpr_str <- indiv_obj$regexpr_str
  
  fit_obj <- readRDS(paste0(model_dir,level,"/",individuals[1],regexpr_str))
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
    fit <- readRDS(paste0(model_dir,level,"/",baboon,regexpr_str))$fit
    clr_ys <- driver::clr(t(fit$Y) + pc)
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
  ggsave(paste0(GP_plot_dir,level,"/all_Sigma_diag_",lrtransform,".png"), plot=p, scale=2, width=16, height=8, units="in", dpi=72)
}

# ====================================================================================================================
# PREDICTION PLOTTING (RIBBONS)
# ====================================================================================================================

# generate ribbon plots for stray::basset fit; assumes predictions have been made in the CLR!!!
#   fit_obj contains the strayfit/bassetfit object and metadata
#   predict object contains the predictions from the strayfit/bassetfit object
#   LR_coord is the logratio coordinate to plot
plot_predictions <- function(fit_obj, predict_obj, LR_coord=1, save_filename=NULL) {
  observations <- fit_obj$X
  clr_ys <- driver::clr(t(fit_obj$Y) + pc)
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
    ylab("LR coord")
  #ylim(c(-10, 10))
  if(is.null(save_filename)) {
    show(p)
  } else {
    ggsave(paste0(GP_plot_dir,level,"/",save_filename,".png"), scale=2, width=12, height=1.5, units="in", dpi=100)
  }
}

plot_ribbons_individuals <- function(individuals, level, timecourse=FALSE, covcor=FALSE, predict_coords=NULL) {
  data <- load_and_filter(level)
  for(baboon in individuals) {
    baboon <<- baboon
    cat(paste0("Visualizing individual '",baboon,"' at level '",level,"'...\n"))
    indiv_data <- subset_samples(data, sname==baboon)
    if(timecourse) {
      cat("\tPlotting timecourse...\n")
      plot_timecourse_phyloseq(indiv_data, paste0(GP_plot_dir,level,"/",baboon,"_timecourse"), gapped=FALSE, 
                               legend=FALSE, legend_level=level)
      plot_timecourse_phyloseq(indiv_data, paste0(GP_plot_dir,level,"/",baboon,"_timecourse"), gapped=TRUE, 
                               legend=FALSE, legend_level=level)
    }
    fit_obj <- readRDS(paste0(model_dir,level,"/",baboon,"_bassetfit.rds"))
    if(covcor) {
      cat(paste0("Plotting element-wise mean covariance/correlation for individual '",baboon,"'...\n"))
      fit.clr <- to_clr(fit_obj$fit)
      Sigma <- fit.clr$Sigma
      meanSigma <- apply(Sigma, c(1,2), mean)
      df <- driver::gather_array(meanSigma, "value", "feature_row", "feature_col")
      p <- ggplot(df, aes(feature_row, feature_col)) +
        geom_tile(aes(fill = value), colour = "white") +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
      ggsave(paste0(GP_plot_dir,level,"/",baboon,"_mean_cov.png"), plot=p, scale=1.5, width=7, height=6, units="in", dpi=72)
      meanSigma_corr <- cov2cor(meanSigma)
      df <- driver::gather_array(meanSigma_corr, "value", "feature_row", "feature_col")
      p <- ggplot(df, aes(feature_row, feature_col)) +
        geom_tile(aes(fill = value), colour = "white") +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
      ggsave(paste0(GP_plot_dir,level,"/",baboon,"_mean_corr.png"), plot=p, scale=1.5, width=7, height=6, units="in", dpi=72)
    }
    if(!is.null(predict_coords) & length(predict_coords) > 0) {
      cat(paste0("Generating predictive plots for individual '",baboon,"'...\n"))
      
      # a dumb hack for now; these need to be global
      se_weight <<- fit_obj$kernelparams$se_weight
      se_sigma <<- fit_obj$kernelparams$se_sigma
      rho_se <<- fit_obj$kernelparams$rho_se
      per_weight <<- fit_obj$kernelparams$per_weight
      per_sigma <<- fit_obj$kernelparams$per_sigma
      rho_per <<- fit_obj$kernelparams$rho_per
      period <<- fit_obj$kernelparams$period
      wn_weight <<- fit_obj$kernelparams$wn_weight
      
      predict_obj <- get_predictions(fit_obj$X, fit_obj$fit, n_samples=100)

      for(coord in predict_coords) {
        plot_predictions(fit_obj, predict_obj, LR_coord=coord, save_filename=paste0(baboon,"_",coord))
      }
    }
  }
}


