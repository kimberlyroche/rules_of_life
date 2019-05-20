library(LaplacesDemon)
library(mvtnorm)
library(matrixsampling)
library(driver)
library(ggplot2)
library(RColorBrewer)

source("include.R")

# ----------------------------------------------------------------------------------------
# (1) simulate some DLM data and fit it (old)
# ----------------------------------------------------------------------------------------

if(FALSE) {
  save_all <- FALSE
  sim <- five_taxa_simulation(T=40, indep_taxa=TRUE, uniform_start=FALSE, noise_scale=1, save_images=FALSE)
  
  cat("Trace true Sigma:",sum(diag(sim$Sigma)),"\n")
  cat("Trace true W:",sum(diag(sim$W)),"\n")
  if(save_all) {
    png(paste0("plots/DLMsim_true_Sigma_",tag,".png"), width=500, height=500)
    image(sim$Sigma)
    dev.off()
  }
  
  censor_vec <- rbinom(nrow(sim$ys), 1, 0.5) # randomly censor to test imputation
  
  fit.f <- fit_filter(sim, censor_vec=censor_vec)
  
  mean.Sigma <- fit.f$Xi/(fit.f$upsilon - ncol(sim$ys) - 1)
  cat("Trace inferred Sigma:",sum(diag(mean.Sigma)),"\n")
  if(save_all) {
    png(paste0("plots/DLMsim_mean_Sigma_",tag,".png"), width=500, height=500)
    image(mean.Sigma)
    dev.off()
  }
  plot_theta_fits(sim, fit_obj.f=fit.f, fit_obj.s=NULL, filename=NULL)
  
  fit.s <- fit_smoother(sim, fit.f, censor_vec=censor_vec)
  
  if(save_all) {
    filename <- paste0("plots/DLMsim_theta_fits_",tag,".png")
  } else {
    filename <- NULL
  }
  plot_theta_fits(sim, fit_obj.f=fit.f, fit_obj.s=fit.s, filename=filename)
  
  # takeaway: can fit a periodic DLM with very few observations if parameters are
  # initialized close to true values and true signal has a very periodic structure
  # the highly structured nature of the signal seems to make it really robust to crappy
  # initial parameter choices; the random walk doesn't have the benefit of that structure!
}

# ----------------------------------------------------------------------------------------
# use the "DUI" and "LEB" (2001) baboon data
# ----------------------------------------------------------------------------------------

sname <- "DUI"

glom_data <- load_glommed_data(level="family", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
# family-agglomerated has replicates; remove these
non_reps <- prune_samples(sample_data(filtered)$sample_status==0, filtered)

D <- 12
indiv_data <- pull_indiv_data(sname, non_reps, subset=D)

# ----------------------------------------------------------------------------------------
# (2a) fit a 1D random walk without gaps; large W here!
# ----------------------------------------------------------------------------------------

data_obj <- list(ys=indiv_data$ys,
                 F=matrix(1, 1, 1),
                 W=diag(1),
                 G=diag(1),
                 upsilon=D+2,
                 Xi=diag(D)*(upsilon-D-1),
                 gamma=1,
                 Sigma=NULL,
                 M.0=matrix(rnorm(D, 1, 1), 1, D),
                 C.0=diag(1)*10)
fit.f <- fit_filter(data_obj, observation_vec=NULL)
fit.s <- fit_smoother(data_obj, fit.f)

filename <- "plots/DLMsim_DUI_1Drandomwalk.png"
filename <- NULL
plot_theta_fits(data_obj, fit.f, fit_obj.s=fit.s, filename=filename) # sanity check

# plot mean Sigma
filename <- "plots/DLMsim_DUI_1Drandomwalk_Sigma.png"
filename <- NULL
plot_mean_Sigma(fit.f, filename=filename)

# ----------------------------------------------------------------------------------------
# (2b) fit a 1D random walk with true between-sample gaps
# ----------------------------------------------------------------------------------------

data_obj <- list(ys=indiv_data$ys,
                 F=matrix(1, 1, 1),
                 W=diag(1)*0.1,
                 G=diag(1),
                 upsilon=D+10,
                 Xi=diag(D)*0.015*(upsilon-D-1),
                 gamma=1,
                 Sigma=NULL,
                 M.0=matrix(rnorm(D, 1, 1), 1, D),
                 C.0=diag(1))
fit.f <- fit_filter(data_obj, observation_vec=indiv_data$observation_vec)
fit.s <- fit_smoother(data_obj, fit.f)

filename <- "plots/DLMsim_DUI_1Drandomwalk_truegaps.png"
filename <- NULL
plot_theta_fits(data_obj, fit_obj.f=NULL, fit_obj.s=fit.s, observation_vec=indiv_data$observation_vec, filename=filename)

filename <- "plots/DLMsim_DUI_1Drandomwalk_gapped_Sigma.png"
filename <- NULL
plot_mean_Sigma(fit.f, filename=filename)

# ----------------------------------------------------------------------------------------
# (3a) fit a FF seasonal model with gapped observations
# use a single harmonic, center Sigma very loosely on HALF the empirical covariance over
# taxa (with the justification that the other half the variance should be learned as)
# system evolution variance and use a discount factor to set W.t; assume gamma = 1
# ----------------------------------------------------------------------------------------

data_obj <- list(ys=indiv_data$ys,
                 F=matrix(c(1, 0), 1, 2),
                 W=NULL,
                 G=build_G(period=365, harmonics=1),
                 upsilon=D+2,
                 Xi=cov(indiv_data$ys)*0.5, # center on empirical covariance
                 gamma=1,
                 Sigma=NULL,
                 M.0=matrix(1, 2, D), # could randomly init too
                 C.0=diag(2)*5)
fit.f <- fit_filter(data_obj, observation_vec=indiv_data$observation_vec, discount=0.98)
fit.s <- fit_smoother(data_obj, fit.f)

filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h1_Theta.png"
filename <- NULL
plot_theta_fits(data_obj, fit_obj.f=NULL, fit_obj.s=fit.s,
                observation_vec=indiv_data$observation_vec, filename=filename)

filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h1_Sigma.png"
filename <- NULL
plot_mean_Sigma(fit.f, filename=filename)

cat("Total variation (empirical):",sum(diag(cov(indiv_data$ys))),"\n")
cat("Total variation (modeled, smoothed):",sum(diag(cov(t(fit.s$etas.t)))),"\n")

# ----------------------------------------------------------------------------------------
# (3b) try a second harmonic
# ----------------------------------------------------------------------------------------

data_obj <- list(ys=indiv_data$ys,
                 F=matrix(c(1, 0, 1, 0), 1, 4),
                 W=NULL,
                 G=build_G(period=365, harmonics=2),
                 upsilon=D+2,
                 Xi=cov(indiv_data$ys)*0.5, # center on empirical covariance
                 gamma=1,
                 Sigma=NULL,
                 M.0=matrix(rnorm(D, 1, 1), 4, D), # could randomly init too
                 C.0=diag(4)*5)
fit.f <- fit_filter(data_obj, observation_vec=indiv_data$observation_vec, discount=0.98)
fit.s <- fit_smoother(data_obj, fit.f)

filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h2_Theta.png"
filename <- NULL
plot_theta_fits(data_obj, fit_obj.f=NULL, fit_obj.s=fit.s,
                observation_vec=indiv_data$observation_vec, filename=filename)

filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h2_Sigma.png"
filename <- NULL
plot_mean_Sigma(fit.f, filename=filename)

cat("Total variation (empirical):",sum(diag(cov(indiv_data$ys))),"\n")
cat("Total variation (modeled, smoothed):",sum(diag(cov(t(fit.s$etas.t)))),"\n")

# ----------------------------------------------------------------------------------------
# (4) do top 10 best sampled
# ----------------------------------------------------------------------------------------

D <- 6
for(sname in best_sampled) {
  cat("Individual",sname,"\n")
  indiv_data <- pull_indiv_data(sname, non_reps, subset_dim=D)
  data_obj <- list(ys=indiv_data$ys,
                   F=matrix(c(1, 0), 1, 2),
                   W=NULL,
                   G=build_G(period=365, harmonics=1),
                   upsilon=D+2,
                   Xi=cov(indiv_data$ys)*0.5, # center on empirical covariance
                   gamma=1,
                   Sigma=NULL,
                   M.0=matrix(1, 2, D), # could randomly init too
                   C.0=diag(2)*5)
  fit.f <- fit_filter(data_obj, observation_vec=indiv_data$observation_vec, discount=0.98)
  fit.s <- fit_smoother(data_obj, fit.f)
  filename <- paste0("plots/DLMsim_",sname,"_1Dperiodic_gaps_h1_Theta.png")
  plot_theta_fits(data_obj, fit_obj.f=NULL, fit_obj.s=fit.s,
                  observation_vec=indiv_data$observation_vec, filename=filename)
  filename <- paste0("plots/DLMsim_",sname,"_1Dperiodic_gaps_h1_Sigma.png")
  plot_mean_Sigma(fit.f, filename=filename)
  cat("Total variation (empirical):",sum(diag(cov(indiv_data$ys))),"\n")
  cat("Total variation (modeled, smoothed):",sum(diag(cov(t(fit.s$etas.t)))),"\n")
}











