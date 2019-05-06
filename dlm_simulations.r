library(LaplacesDemon)
library(mvtnorm)
library(matrixsampling)
library(driver)
library(ggplot2)
library(RColorBrewer)

source("include.R")

save_all <- FALSE

# ----------------------------------------------------------------------------------------
# (1) simulate some DLM data and fit it
# ----------------------------------------------------------------------------------------

sim <- five_taxa_simulation(T=40, indep_taxa=TRUE, uniform_start=FALSE, save_images=FALSE)

cat("Trace true Sigma:",sum(diag(sim$Sigma)),"\n")
if(save_all) {
  png(paste0("plots/DLMsim_true_Sigma_",tag,".png"), width=500, height=500)
  image(sim$Sigma)
  dev.off()
}

censor_vec <- rbinom(nrow(sim$ys), 1, 0.5) # randomly censor half to test imputation

fit.f <- fit_filter(sim, censor_vec=censor_vec)

mean.Sigma <- fit.f$Xi/(fit.f$upsilon - ncol(sim$ys) - 1)
cat("Trace inferred Sigma:",sum(diag(mean.Sigma)),"\n")
if(save_all) {
  png(paste0("plots/DLMsim_mean_Sigma_",tag,".png"), width=500, height=500)
  image(mean.Sigma)
  dev.off()
}

fit.s <- fit_smoother(sim, fit.f, censor_vec=censor_vec)

filename <- paste0("plots/DLMsim_theta_fits_",tag,".png")
plot_theta_fits(sim, fit_obj.f=fit.f, fit_obj.s=fit.s, filename=NULL)

# ----------------------------------------------------------------------------------------
# use the "DUI" (2001) baboon data
# ----------------------------------------------------------------------------------------

glom_data <- load_glommed_data(level="family", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
# family-agglomerated has replicates; remove these
non_reps <- prune_samples(sample_data(filtered)$sample_status==0, filtered)
pruned <- prune_samples(sample_data(non_reps)$sname=="DUI", non_reps)
pruned <- prune_samples((sample_data(pruned)$collection_date > "2001-10-01") &
                          (sample_data(pruned)$collection_date < "2002-11-30"), pruned)
cat("Real data set has",nsamples(pruned),"samples and",ntaxa(pruned),"taxa\n")

counts <- otu_table(pruned)@.Data + 0.65 # samples (rows) x taxa (columns)
ys <- driver::clr(counts)

# ----------------------------------------------------------------------------------------
# (2a) fit a 1D random walk without gaps
# ----------------------------------------------------------------------------------------

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)
G <- diag(1)
upsilon <- D + 2
Xi <- diag(D)*(upsilon-D-1)
gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=NULL)
fit.s <- fit_smoother(data_obj, fit.f, censor_vec=NULL)

plot_theta_fits(data_obj, fit.f, fit_obj.s=fit.s, filename="plots/DLMsim_DUI_1Drandomwalk.png") # sanity check

# plot mean Sigma
mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
if(save_all) {
  png("plots/DLMsim_DUI_1Drandomwalk_Sigma.png")
  image(mean.Sigma)
  dev.off()
}

# ----------------------------------------------------------------------------------------
# (2b) fit a 1D random walk with a big gap in the middle
# ----------------------------------------------------------------------------------------

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)
G <- diag(1)
upsilon <- D + 2
Xi <- diag(D)*(upsilon-D-1)
gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W

observation_vec <- c(1:15, 26:40) # 10-observation gap in the middle

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
fit.s <- fit_smoother(data_obj, fit.f, censor_vec=NULL)

plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, filename="plots/DLMsim_DUI_1Drandomwalk_middlegap.png")
# check the means, not the Theta.t draws
#fit_obj$Thetas.t <- fit_obj$Ms.t
#plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=NULL, filename=NULL)

# ----------------------------------------------------------------------------------------
# (2c) fit a 1D random walk with true between-sample gaps
# ----------------------------------------------------------------------------------------

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)
G <- diag(1)
upsilon <- D + 2
Xi <- diag(D)*(upsilon-D-1)
gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W

dates <- sample_data(pruned)$collection_date
min_d <- dates[1]
observation_vec <- sapply(dates, function(x) { round(difftime(x, min_d, units="days"))+1 } )

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
fit.s <- fit_smoother(data_obj, fit.f, censor_vec=NULL)

plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, filename="plots/DLMsim_DUI_1Drandomwalk_truegaps.png")
# check the means, not the Theta.t draws
#fit_obj$Thetas.t <- fit_obj$Ms.t
#plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=NULL, filename="plots/DLMsim_DUI_1Drandomwalk_truegaps_mean.png")

# plot mean Sigma
mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
if(save_all) {
  png("plots/DLMsim_DUI_1Drandomwalk_gapped_Sigma.png")
  image(mean.Sigma)
  dev.off()
}

# ----------------------------------------------------------------------------------------
# (...) fit a random walk model with a 2-component state
# ----------------------------------------------------------------------------------------

# TO-DO

# ----------------------------------------------------------------------------------------
# (...) fit a FF seasonal model with gapped observations
# ----------------------------------------------------------------------------------------

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)*0.1
G <- build_G(period=365)
upsilon <- D + 2
Xi <- diag(D)*(upsilon-D-1)*0.1
gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W

dates <- sample_data(pruned)$collection_date
min_d <- dates[1]
observation_vec <- sapply(dates, function(x) { round(difftime(x, min_d, units="days"))+1 } )

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec)

plot_theta_fits(data_obj, fit.f, fit_obj.s=NULL, filename=NULL) # sanity check





