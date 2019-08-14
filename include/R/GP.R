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
# not currently in use
WHITENOISE <- function(X, sigma=1, jitter=1e-10) {
  dist <- as.matrix(dist(t(X)))
  G <- diag(ncol(dist))*sigma^2 + jitter*diag(ncol(dist))
  return(G)
}

# ====================================================================================================================
# STRAY WRAPPERS
# ====================================================================================================================

# fit Gaussian process to a single baboon series using stray::basset
fit_GP <- function(baboon, level, se_weight=2, per_weight=0.25, dd_se=90, save_append="",
                   date_lower_limit=NULL, date_upper_limit=NULL, alr_ref=NULL) {
  cat(paste0("Fitting stray::basset with with parameters:\n",
             "\tbaboon=",baboon,"\n",
             "\tlevel=",level,"\n",
             "\tSE kernel weight=",se_weight,"\n",
             "\tPER kernel weight=",per_weight,"\n",
             "\tdd_se=",dd_se,"\n"))

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
    (1e-8)*diag(ncol(X)) # pretty arbitrary
  
  D <- nrow(Y)
  N <- ncol(Y)
  
  # stray uses the D^th element as the ALR reference by default
  # do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }
  
  # ALR prior covariance
  upsilon <- D-1+10 # lesser certainty
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  # take diag as covariance over log abundances; scale 1 is probably shit here but
  Xi <- GG%*%(diag(D)*1)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
  
  alr_ys <- driver::alr((t(Y)+0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
  fit <- stray::basset(Y, observations, upsilon, Theta, Gamma, Xi)
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
  
  # chop down for size savings since /data/mukherjeelab is full as shit
  # can't seem to pass desired sample number to stray with any effect; debug this eventually
  fit_obj$fit$iter <- 100
  fit_obj$fit$Eta <- fit_obj$fit$Eta[,,1:fit_obj$fit$iter]
  fit_obj$fit$Lambda <- fit_obj$fit$Lambda[,,1:fit_obj$fit$iter]
  fit_obj$fit$Sigma <- fit_obj$fit$Sigma[,,1:fit_obj$fit$iter]
  
  saveRDS(fit_obj, paste0(model_dir,level,"/",baboon,"_bassetfit",save_append,".rds"))
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

# embed posterior samples of (which_measure) using MDS and the appropriate distance metric
#   which_measure: Sigma | Lambda
embed_posteriors <- function(level, which_measure="Sigma") {
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
              mean_x1=mean(x3), mean_x2=mean(x4),
              mean_x1=mean(x5), mean_x2=mean(x6))
  saveRDS(df_centroids, paste0(GP_plot_dir,level,"/",which_measure,"_ordination_centroids.rds"))
}

# ====================================================================================================================
# EMBEDDING VISUALIZATION FNS.
# ====================================================================================================================

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
  # axes 1 & 2
  plot_axes(df, df_centroids, "x1", "x2", label_type, legend=legend)
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

# plot on- or off-diagonal elements from covariance matrices for most extreme individuals in
# chosen coordinate of BASELINE embedding; this is a terrible visualization
#   coordinate is 1...K
#   level is taxonomic level (e.g. family)
#   no_indiv are # individuals at each extreme to visualize
plot_extreme_Sigma <- function(coordinate, level, no_indiv=10, save_filename="extreme_Sigma") {
  df_centroids <- readRDS(paste0(GP_plot_dir,level,"/Sigma_ordination_centroids.rds"))
  min_sort <- df_centroids %>% arrange(get(coordinate))
  min_cohort <- as.vector(unlist(min_sort[1:no_indiv,"labels"]))
  max_sort <- df_centroids %>% arrange(desc(get(coordinate)))
  max_cohort <- as.vector(unlist(max_sort[1:no_indiv,"labels"]))
  
  df_on <- data.frame(sample=c(), feature=c(), value=c())
  df_off <- data.frame(sample=c(), feature=c(), value=c())
  for(baboon in c(min_cohort, max_cohort)) {
    cat("Loading individual",baboon,"\n")
    Sigma <- readRDS(paste0(model_dir,level,"/",baboon,"_bassetfit.rds"))$fit$Sigma
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
  ggsave(paste0(GP_plot_dir,level,"/",save_filename,"_ondiag.png"), plot=p, scale=2, width=6, height=6, units="in", dpi=72)
  p <- ggplot(df_off, aes(feature, sample)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred")
  ggsave(paste0(GP_plot_dir,level,"/",save_filename,"_offdiag.png"), plot=p, scale=2, width=10, height=6, units="in", dpi=72)
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
  ggsave(paste0(GP_plot_dir,level,"/all_Sigma_diag_",lrtransform,".png"), plot=p, scale=2, width=16, height=8, units="in", dpi=72)
}


