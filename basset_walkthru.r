library(stray)
library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)

# X is Q x N as in other kernels
# bandwidth: rho as chosen gives antiphase observations a correlation of ~0.1
PER <- function(X, sigma=1, rho=1, period = 24, jitter = 1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

# returns fit (CLR)
fit_to_reference <- function(Gamma){
  data(mallard) # loads "ps" phyloseq obj
  
  # massage data
  ps <- prune_samples(sample_data(ps)$Vessel==1, ps) # vessel 1
  ps <- prune_samples((sample_data(ps)$time > "2015-11-20 15:00:00 UTC") &
                        (sample_data(ps)$time < "2015-11-25 16:00:00 UTC"), ps) # hourly samples only
  o <- order(sample_data(ps)$time) # get the order of samples
  otu_table(ps) <- otu_table(ps)[o,] # extract counts in this order (samples are rows!)
  sample_data(ps) <- sample_data(ps)[o,] # extract sample data in this order (samples are rows!)
  Y <- t(as(otu_table(ps), "matrix")) # extract counts; D x N
  rownames(Y) <- taxa_names(ps) # format: seq_*

  D <- ntaxa(ps)
  N <- nrow(sample_data(ps))
  
  upsilon <- D-1+3
  Xi <- matrix(.4, D-1, D-1)
  diag(Xi) <- 1
  Theta <- function(X) matrix(0, D-1, ncol(X))
  
  X <- as.numeric(sample_data(ps)$time) # hour identifier
  X <- t((X-min(X)) / 3600) # readable hours 0..115
  
  fit <- stray::basset(Y, X, upsilon, Theta, Gamma, Xi)
  fit.clr <- to_clr(fit)

  family_names <- as(tax_table(ps)[,"Family"], "vector")

  return(list(Y=Y, X=X, fit=fit.clr, labels=as(tax_table(ps)[,"Family"], "vector")))
}

fit_to_baboon <- function(Gamma, sname, start_date_str, end_date_str) {
  source("include.R")
  glom_data <- load_glommed_data(level="family", replicates=TRUE)
  filtered <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
  # family-agglomerated has replicates; remove these
  non_reps <- prune_samples(sample_data(filtered)$sample_status==0, filtered)
  cat("Filtered down to",ntaxa(non_reps),"taxa\n")
  
  pruned <- prune_samples(sample_data(non_reps)$sname==sname, non_reps)
  pruned <- prune_samples((sample_data(pruned)$collection_date > start_date_str) &
                        (sample_data(pruned)$collection_date < end_date_str), pruned)
  cat("Using",phyloseq::nsamples(pruned),"samples from",sname,"\n")
  tax <- tax_table(pruned)
  o <- order(sample_data(pruned)$collection_date)
  otu_table(pruned) <- otu_table(pruned)[o,]
  sample_data(pruned) <- sample_data(pruned)[o,]
  Y <- t(as(otu_table(pruned), "matrix")) # dimensions D x N
  readable_rownames <- c()
  base_names_used <- list()
  for(t in taxa_names(pruned)) {
    # we need these identifiers to be unique, so if family is unknown (frequently)
    # backtrack to find some usable identifier
    tax_lvl <- 5
    tax_label <- as(tax[t,], "vector")[tax_lvl]
    # if name is already taken, just append a counter to it
    if(tax_label %in% names(base_names_used)) {
      base_names_used[[tax_label]] = base_names_used[[tax_label]] + 1
      tax_label <- paste0(tax_label, " (", base_names_used[[tax_label]], ")")
    } else {
      # otherwise if it's NA, move up phylogenetic tree until we hit something concrete
      while(is.na(tax_label)) {
        tax_lvl <- tax_lvl - 1
        tax_label <- as(tax[t,], "vector")[tax_lvl]
      }
      if(tax_lvl == 0 && is.na(tax_label)) {
        tax_label <- "unknown"
      }
      if(tax_label %in% names(base_names_used)) {
        base_names_used[[tax_label]] = base_names_used[[tax_label]] + 1
        tax_label <- paste0(tax_label, " (", base_names_used[[tax_label]], ")")
      } else {
        base_names_used[[tax_label]] = 1
      }
    }
    readable_rownames <- c(readable_rownames, tax_label)
  }
  rownames(Y) <- readable_rownames
  
  D <- ntaxa(pruned)
  N <- nrow(sample_data(pruned))
  
  upsilon <- D-1+3
  Xi <- matrix(.4, D-1, D-1)
  diag(Xi) <- 1
  Theta <- function(X) matrix(0, D-1, ncol(X))
  
  X_raw <- sample_data(pruned)$collection_date
  # sorted, so min date is first
  X <- sapply(X_raw, function(Z) { round(difftime(Z, X_raw[1], units="days")) } )
  dim(X) <- c(1, length(X))
  
  fit <- stray::basset(Y, X, upsilon, Theta, Gamma, Xi)
  fit.clr <- to_clr(fit)
  return(list(Y=Y, X=X, fit=fit.clr, labels=readable_rownames))
}

get_predictions <- function(X, fit.clr){
  cat("Predicting from 1 to",max(X),"\n")
  X_predict <- t(1:(max(X))) # time point range, fill in any missing
  predicted <- predict(fit.clr, X_predict, jitter=1) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

plot_fit <- function(Y, X, labels, Y_predict, X_predict, save_filename){
  Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>% 
    gather_array(mean, coord, sample) %>% 
    mutate(time = X[1,sample], 
           coord = paste0("CLR(", labels[coord],")"))
  
  predicted_tidy <- gather_array(Y_predict, val, coord, sample, iter) %>% 
    mutate(time = X_predict[1,sample]) %>% # just adds time point ID
    filter(!is.na(val)) %>% # kill NAs
    group_by(time, coord) %>% # there are D x N combos of time and coordinate (log ratio)
    summarise_posterior(val, na.rm=TRUE) %>% # gets standard quantiles for each coord x timepoint combo
    ungroup() %>% # just undoes this group x coord binding
    mutate(coord = paste0("CLR(", labels[coord],")")) # replace coord with readable family name
  
  p <- ggplot(predicted_tidy, aes(x = time, y=mean)) +
    geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(color="blue") +
    geom_point(data = Y_clr_tidy, alpha=0.5) +
    facet_wrap(~coord, scales="free_y") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle=45))
  ggsave(paste0("plots/",save_filename,".png",sep=""), plot=p, scale=1.5, width=8, height=5, units="in")
}

if(FALSE) {
  # fit to MALLARD data
  Gamma <- function(X) SE(X) # squared exponential only
  fit.obj <- fit_to_reference(Gamma)
  predict.obj <- get_predictions(fit.obj$X, fit.obj$fit)
  plot_fit(fit.obj$Y, fit.obj$X, fit.obj$labels, predict.obj$Y_predict, predict.obj$X_predict)
}

# baboon
#Gamma <- function(X) PER(X, period=365) # periodic only
Gamma <- function(X) SE(X) # squared exponential only
for(indiv in c("DUI", "LEB")) {
  fit.obj <- fit_to_baboon(Gamma, indiv, "2001-10-01", "2002-11-30") # wet season through a dry season
  predict.obj <- get_predictions(fit.obj$X, fit.obj$fit)
  plot_fit(fit.obj$Y, fit.obj$X, fit.obj$labels, predict.obj$Y, predict.obj$X, paste0("basset_SE_",indiv))
}










