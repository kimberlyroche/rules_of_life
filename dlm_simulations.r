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

saved_ys <- driver::clr(counts)
ys <- saved_ys
ys <- ys[,1:10]

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
fit.s <- fit_smoother(data_obj, fit.f)

if(save_all) {
  filename <- "plots/DLMsim_DUI_1Drandomwalk.png"
} else {
  filename <- NULL
}
plot_theta_fits(data_obj, fit.f, fit_obj.s=fit.s, filename=filename) # sanity check

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

saved_ys <- driver::clr(counts)
ys <- saved_ys
ys <- ys[,1:10]

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)
G <- diag(1)
upsilon <- D + 2
Xi <- diag(D)*(upsilon-D-1)
gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W

observation_vec <- c(1:10, 21:30) # 10-observation gap in the middle

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
fit.s <- fit_smoother(data_obj, fit.f, censor_vec=NULL)

if(save_all) {
  filename <- "plots/DLMsim_DUI_1Drandomwalk_middlegap.png"
} else {
  filename <- NULL
}
plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, filename=filename)

# check the means, not the Theta.t draws
fit.f$Thetas.t <- fit.f$Ms.t
fit.s$Thetas.t <- fit.s$Ms.t
plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, filename=NULL)

mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
cat("Trace inferred Sigma:",sum(diag(mean.Sigma)),"\n")
image(mean.Sigma)

# takeaway: the means track pretty well where they have feedback and more or less
# flatline where they don't; with a gap in available observations in the middle, 
# smoothing definitely improves on the filtered estimate, which is SHOULD; small values
# for W give lesser/slower adaptation and smaller Sigma but really don't change structure
# of Sigma

# there's a tradeoff -- no missing: large W make for very tight fits, small W make for
# very smooth fits; gaps: small W keep the thetas from wandering too much in the gaps
# (in either case the mean looks good)

# ----------------------------------------------------------------------------------------
# (2c) fit a 1D random walk with true between-sample gaps
# ----------------------------------------------------------------------------------------

saved_ys <- driver::clr(counts)
ys <- saved_ys
ys <- ys[,1:4]

D <- ncol(ys)

F <- matrix(1, 1, 1)
W <- diag(1)*0.1
G <- diag(1)
upsilon <- D + 10
Xi <- diag(D)*0.025*(upsilon-D-1) # setting these parameters in such a way that the taxa
                                # damn near totally independent and the system and
                                # observation noise are very small seems to work well

cat("Relative scale of W, Sigma:\n")
cat("Tr W:",sum(diag(W)),"\n")
Sigma_prior <- Xi/(upsilon - D - 1)
cat("Tr Sigma:",sum(diag(Sigma_prior)),"\n")

gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 1, D)
C.0 <- W*10

dates <- sample_data(pruned)$collection_date
min_d <- dates[1]
observation_vec <- sapply(dates, function(x) { round(difftime(x, min_d, units="days"))+1 } )

data_obj <- list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
fit.s <- fit_smoother(data_obj, fit.f)

if(save_all) {
  filename <- "plots/DLMsim_DUI_1Drandomwalk_truegaps.png"
} else {
  filename <- NULL
}
plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, observation_vec=observation_vec, filename=filename)

# check the means, not the Theta.t draws
#alt_fit.f <- fit.f
#alt_fit.s <- fit.s
#alt_fit.f$Thetas.t <- alt_fit.f$Ms.t
#alt_fit.s$Thetas.t <- alt_fit.s$Ms.t
#plot_theta_fits(data_obj, fit_obj.f=alt_fit.f, fit_obj.s=alt_fit.s, observation_vec=observation_vec, filename=NULL)

# plot mean Sigma
mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
if(save_all) {
  png("plots/DLMsim_DUI_1Drandomwalk_gapped_Sigma.png")
  image(mean.Sigma)
  dev.off()
}

cat("Total variation (empirical):",sum(diag(cov(ys))),"\n")
cat("Total variation (modeled, filtered):  ",(sum(diag(cov(t(fit.f$Thetas.t[1,,])))) + gamma*sum(diag(mean.Sigma))),"\n")
cat("Total variation (modeled, smoothed):  ",(sum(diag(cov(t(fit.s$Thetas.t[1,,])))) + gamma*sum(diag(mean.Sigma))),"\n")

# ----------------------------------------------------------------------------------------
# (3a) fit a FF seasonal model with gapped observations
# use a single harmonic, center Sigma very loosely on HALF the empirical covariance over
# taxa (with the justification that the other half the variance should be learned as)
# system evolution variance and use a discount factor to set W.t; assume gamma = 1
# ----------------------------------------------------------------------------------------

saved_ys <- driver::clr(counts)
ys <- saved_ys
ys <- ys[,1:6]

D <- ncol(ys)

F <- matrix(c(1, 0), 1, 2)
G <- build_G(period=365, harmonics=1)
upsilon <- D + 2 # very weak centering
Xi <- cov(ys)*0.5 # center on empirical covariance

cat("Relative scale of Sigma (trace):",sum(diag(Xi)),"\n")

gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 2, D)
C.0 <- diag(2)

dates <- sample_data(pruned)$collection_date
min_d <- dates[1]
observation_vec <- sapply(dates, function(x) { round(difftime(x, min_d, units="days"))+1 } )

data_obj <- list(ys=ys, F=F, W=NULL, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec, discount=0.98)
fit.s <- fit_smoother(data_obj, fit.f)

if(save_all) {
  filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h1_Theta.png"
} else {
  filename <- NULL
}
plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, observation_vec=observation_vec, filename=filename)

# plot mean Sigma
mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
if(save_all) {
  png("plots/DLMsim_DUI_1Dperiodic_gaps_h1_Sigma.png")
  image(mean.Sigma)
  dev.off()
}

cat("Total variation (empirical):",sum(diag(cov(ys))),"\n")
T <- dim(fit.s$Thetas.t)[3]
etas.f <- matrix(0, D, T)
etas.s <- matrix(0, D, T)
for(i in 1:T) {
  etas.f[,t] <- F%*%fit.f$Thetas.t[,,t]
  etas.s[,t] <- F%*%fit.s$Thetas.t[,,t]
}
cat("Total variation (modeled, filtered):  ",(sum(diag(cov(etas.f))) + gamma*sum(diag(mean.Sigma))),"\n")
cat("Total variation (modeled, smoothed):  ",(sum(diag(cov(etas.s))) + gamma*sum(diag(mean.Sigma))),"\n")

# ----------------------------------------------------------------------------------------
# (3b) periodic + half harmonic
# ----------------------------------------------------------------------------------------

saved_ys <- driver::clr(counts)
ys <- saved_ys
ys <- ys[,1:6]

D <- ncol(ys)

F <- matrix(c(1, 0, 1, 0), 1, 4)
G <- build_G(period=365, harmonics=2)
upsilon <- D + 2 # very weak centering
Xi <- cov(ys)*0.5 # center on empirical covariance

cat("Relative scale of Sigma (trace):",sum(diag(Xi)),"\n")

gamma <- 1
M.0 <- matrix(rnorm(D, 1, 1), 4, D)
C.0 <- diag(4)

dates <- sample_data(pruned)$collection_date
min_d <- dates[1]
observation_vec <- sapply(dates, function(x) { round(difftime(x, min_d, units="days"))+1 } )

data_obj <- list(ys=ys, F=F, W=NULL, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=NULL, M.0=M.0, C.0=C.0)
fit.f <- fit_filter(data_obj, observation_vec=observation_vec, discount=0.98)
fit.s <- fit_smoother(data_obj, fit.f)

if(save_all) {
  filename <- "plots/DLMsim_DUI_1Dperiodic_gaps_h2_Theta.png"
} else {
  filename <- NULL
}
plot_theta_fits(data_obj, fit_obj.f=fit.f, fit_obj.s=fit.s, observation_vec=observation_vec, filename=filename)

# plot mean Sigma
mean.Sigma <- fit.f$Xi/(fit.f$upsilon - D - 1)
if(save_all) {
  png("plots/DLMsim_DUI_1Dperiodic_gaps_h2_Sigma.png")
  image(mean.Sigma)
  dev.off()
}

cat("Total variation (empirical):",sum(diag(cov(ys))),"\n")
T <- dim(fit.s$Thetas.t)[3]
etas.f <- matrix(0, D, T)
etas.s <- matrix(0, D, T)
for(i in 1:T) {
  etas.f[,t] <- F%*%fit.f$Thetas.t[,,t]
  etas.s[,t] <- F%*%fit.s$Thetas.t[,,t]
}
cat("Total variation (modeled, filtered):  ",(sum(diag(cov(etas.f))) + gamma*sum(diag(mean.Sigma))),"\n")
cat("Total variation (modeled, smoothed):  ",(sum(diag(cov(etas.s))) + gamma*sum(diag(mean.Sigma))),"\n")

