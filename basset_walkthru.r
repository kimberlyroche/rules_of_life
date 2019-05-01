library(stray)
library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)

# X is Q x N as in other kernels
# bandwidth chosen as in SE
PER <- function(X, sigma = 1, rho = median(as.matrix(dist(t(X)))), period = 24, jitter = 1e-10){
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
  return(list(Y=Y, X=X, fit=fit.clr, labels=as(tax_table(ps)[,"Family"], "vector")))
}

fit_to_ACA <- function(Gamma) {
  source("include.R")
  glom_data <- load_glommed_data(level="genus", replicates=FALSE) # to do: smash this to family
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
  
  indiv <- "ACA"
  pruned <- prune_samples(sample_data(filtered)$sname==indiv, filtered) # ACA
  pruned <- prune_samples((sample_data(pruned)$collection_date > "2012-01-01") &
                        (sample_data(pruned)$collection_date < "2019-01-01"), pruned) # hourly samples only
  tax <- tax_table(pruned)
  o <- order(sample_data(pruned)$collection_date)
  otu_table(pruned) <- otu_table(pruned)[o,]
  sample_data(pruned) <- sample_data(pruned)[o,]
  Y <- t(as(otu_table(pruned), "matrix")) # dimensions D x N
  readable_rownames <- c()
  for(t in taxa_names(pruned)) {
    # use family/genus
    readable_rownames <- c(readable_rownames, paste(as(tax[t,], "vector")[5:6],collapse="/"))
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

get_predictions <- function(Y, X, fit.clr){
  X_predict <- t(1:(max(X))) # time point range, fill in any missing
  predicted <- predict(fit.clr, X_predict, jitter=1) # predicts samples from the posterior (default = 2000)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

plot_fit <- function(Y, X, labels, Y_predict, X_predict){
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
  p
}

Gamma <- function(X) SE(X, sigma=5, rho=10) # squared exponential only
#Gamma <- function(X) PER(X, sigma=1, rho=5, period=24) # periodic only
#Gamma <- function(X) PER(X, sigma=1, rho=5, period=24) + diag(ncol(X))*0.01 + LINEAR(X, sigma=1) # Create partial function 

#fit.obj <- fit_to_reference(Gamma) # priors specified as function
fit.obj <- fit_to_ACA(Gamma)
predict.obj <- get_predictions(fit.obj$Y, fit.obj$X, fit.obj$fit)
plot_fit(fit.obj$Y[1:5,], fit.obj$X, fit.obj$labels[1:5], predict.obj$Y[1:5,,], predict.obj$X)











