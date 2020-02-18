rm(list=ls())

library(matrixsampling)
library(stray)
library(Rcpp)
library(driver)
library(ggplot2)
library(gridExtra)

sourceCpp("include/cpp/Riemann_dist.cpp")

# question of interest: do similar compositions always give similar dynamics?
#                       do different compositions always give different dynamics?

# in all cases, simulate data (up to counts), fit model, and assess Sigma true vs.
# Sigma posteriors via PCoA

# setup 1: similar base compositions (Theta), different dynamics (Sigma)

# setup 2: different base compositions (Theta), same dynamics (Sigma)

setup <- 1

D <- 20
N <- 20
K <- 10 # number of individuals to simulate 
n_samples <- 100 # posterior samples

# time index
X <- matrix(1:N, 1, N)

# same Theta (function) and Gamma
Theta_fn <- function(X) {
  matrix(Thetas[,X[2,1],drop=F], D-1, ncol(X))
}

Gamma_fn <- function(X) SE(X[1,,drop=F], sigma=1, rho=5)

# same hyperparameters
upsilon <- D-1+3 # lesser certainty
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

# potentially simulate K "baselines"
Thetas <- matrix(NA, D-1, K)
# potentially simulate K "dynamics", sets of counts, etc.
Sigmas <- matrix(NA, D-1, K*(D-1)*(n_samples+1)) # true + posterior sampled
# matrix `Sigmas` this will have the form ((D-1) x (D-1)*(n_samples+1)*K), where below M = n_samples:
#
# ----------------------------------------------------------------------------------------------------------- -
# |    indiv 1    |    indiv 1    |    indiv 1    |               |    indiv 1    |    indiv 2    |
# |     truth     | fit sample 1  | fit sample 2  |      ...      | fit sample M  |     truth     |      ...
# | (D-1) x (D-1) | (D-1) x (D-1) | (D-1) x (D-1) |               | (D-1) x (D-1) | (D-1) x (D-1) |
# ----------------------------------------------------------------------------------------------------------- -
#
# Sigma_list <- list() # this is just redundant (temporary) storage of true Sigmas so I can double check my indexing

# convenience function
get_Sigma_samples <- function(k, all_samples, D, n_samples) {
  samples <- NULL
  truth_idx1 <- (k-1)*(D-1)*(n_samples+1) + 1 + (D-1) - 1
  fit_idx1 <- truth_idx1 + 1
  fit_idx2 <- fit_idx1 + (D-1) - 1
  for(j in 1:n_samples) {
    Sigma_sample <- Sigmas[,fit_idx1:fit_idx2]
    if(is.null(samples)) {
      samples <- Sigma_sample
    } else {
      samples <- cbind(samples, Sigma_sample)
    }
    fit_idx1 <- fit_idx2 + 1
    fit_idx2 <- fit_idx1 + (D-1) - 1
  }
  return(samples)
}

labels <- data.frame(status=c(), k=c())

# construct `Thetas` in a first pass, since Theta_fn will need to reference it
for(k in 1:K) {
  if(setup == 1) {
    # all individuals will have the same baseline composition but different (random) dynamics
    if(k == 1) {
      Theta <- rnorm(D-1, 0, 1.5)
    } else {
      Theta <- Thetas[,1]
    }
    Thetas[,k] <- Theta
  } else {
    # all individuals will have the different (random) baseline compositions but the same dynamics
    Theta <- rnorm(D-1, 0, 2)
    Thetas[,k] <- Theta
  }
}

fits <- list()
fits_MAP <- list()
Etas <- list()
Lambdas <- list()
Xs <- list()
for(k in 1:K) {
  cat("Generating data set",k,"...\n")
  if(setup == 1) {
    Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  } else {
    if(k == 1) {
      Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
    } else {
      Sigma <- Sigmas[,1:(D-1)]
    }
  }

  # terrible indexing math! but you only have to figure it out once...
  truth_idx1 <- (k-1)*(D-1)*(n_samples+1)+1
  truth_idx2 <- truth_idx1 + (D-1) - 1
  Sigmas[,truth_idx1:truth_idx2] <- Sigma
  #Sigma_list[[k]] <- Sigma
  Xs[[k]] <- rbind(X, k)
  Lambda <- rmatrixnormal(1, Theta_fn(Xs[[k]]), Sigma, Gamma_fn(rbind(X, k)))[,,1]
  Lambdas[[k]] <- Lambda
  eta <- rmatrixnormal(1, Lambda, Sigma, diag(N))[,,1] # D-1 x N
  Etas[[k]] <- eta
  pi <- matrix(NA, D, N)
  Y <- matrix(NA, D, N)
  for(j in 1:N) {
    pi[,j] <- driver::alrInv(eta[,j])
    Y[,j] <- rmultinom(1, rpois(1, 10000), pi[,j])
  }
  
  # sanity check; plot as log relative abundance paths
  df <- gather_array(eta, "value", "taxon", "sample")
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df, aes(x=sample, y=value, color=taxon)) +
    geom_path() +
    theme(legend.position = "none")
  if(setup == 1) {
    ggsave(paste0("output/same_baseline_diff_dynamics_path_",k,".png"), p, units="in", dpi=150, height=4, width=10)
  } else {
    ggsave(paste0("output/diff_baseline_same_dynamics_path_",k,".png"), p, units="in", dpi=150, height=4, width=10)
  }
  
  # sanity check; plot as proportions
  df <- gather_array(pi, "proportion", "taxon", "sample")
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df, aes(x=sample, y=proportion, fill=taxon)) + 
    geom_bar(position="fill", stat="identity") +
    theme(legend.position = "none")
  if(setup == 1) {
    ggsave(paste0("output/same_baseline_diff_dynamics_barplot_",k,".png"), p, units="in", dpi=150, height=4, width=10)
  } else {
    ggsave(paste0("output/diff_baseline_same_dynamics_barplot_",k,".png"), p, units="in", dpi=150, height=4, width=10)
  }

  fit <- basset(Y, rbind(X, k), upsilon, Theta_fn, Gamma_fn, Xi, n_samples=0, ret_mean=TRUE)
  fits_MAP[[k]] <- fit
  fit <- basset(Y, rbind(X, k), upsilon, Theta_fn, Gamma_fn, Xi, n_samples=n_samples)
  fits[[k]] <- fit
  fit_idx1 <- truth_idx2 + 1
  fit_idx2 <- fit_idx1 + (D-1) - 1
  for(j in 1:n_samples) {
    Sigmas[,fit_idx1:fit_idx2] <- fit$Sigma[,,j]
    fit_idx1 <- fit_idx2 + 1
    fit_idx2 <- fit_idx1 + (D-1) - 1
  }

  labels <- rbind(labels, data.frame(status=c("truth", rep("fitted", n_samples)), k=k))
}

# let's spike in some noise around the baseline and see where that ends up on the ordination
# to help us understand the scale of the distances here
do_spike_in <- FALSE

if(do_spike_in) {
  spike_in <- matrix(NA, D-1, (D-1)*(n_samples+1))
  spike_in[,1:(D-1)] <- Sigmas[,1:(D-1)]
  for(i in 1:n_samples) {
    offset <- (D-1)*i
    left <- offset + 1
    right <- offset + (D-1)
    upsilon <- (D-1) + 2 + 50
    spike_in[,left:right] <- rinvwishart(1, upsilon, spike_in[,1:(D-1)]*(upsilon - (D-1) - 1))[,,1]
  }
  Sigmas <- cbind(Sigmas, spike_in)
  labels <- rbind(labels, data.frame(status=c("truth", rep("fitted", n_samples)), k=rep(0, n_samples+1)))

  plot(density(spike_in), ylim=c(0,3))
  for(k in 1:K) {
    lines(density(get_Sigma_samples(k, all_samples, D, n_samples)), col="red")
  }
}

# not fast
# for D=20 / N=20 / K=10 / n_samples=100 this takes about 1 min.
n_Sigma <- ncol(Sigmas)/(D-1)
d <- matrix(NA, n_Sigma, n_Sigma)
for(i in 1:n_Sigma) {
  for(j in 1:n_Sigma) {
    if(i <= j) {
      d[i,j] <- Riemann_dist_pair(Sigmas[,((i-1)*(D-1)+1):(i*(D-1))], Sigmas[,((j-1)*(D-1)+1):(j*(D-1))])
    } else {
      d[i,j] <- d[j,i]
    }
  }
}

# embed this inferred Sigmas and the true Sigmas to see if the inferred stuff is "too similar"
embedding <- cmdscale(d, k=4)
df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                 x3=embedding[,3], x4=embedding[,4],
                 type=as.factor(labels$status), indiv=as.factor(labels$k))
fitted_point_sz <- 2
truth_point_size <- 4
p1 <- ggplot() +
  geom_point(data=df[df$type=="fitted",], aes(x=x1, y=x2, color=indiv), size=fitted_point_sz) +
  geom_point(data=df[df$type=="truth",], aes(x=x1, y=x2), color="black", size=(truth_point_size*1.3), shape=17) +
  geom_point(data=df[df$type=="truth",], aes(x=x1, y=x2, color=indiv), size=truth_point_size, shape=17) +
  xlab("PCoA 1") +
  ylab("PCoA 2")

p1

p2 <- ggplot() +
  geom_point(data=df[df$type=="fitted",], aes(x=x3, y=x4, color=indiv), size=fitted_point_sz) +
  geom_point(data=df[df$type=="truth",], aes(x=x3, y=x4), color="black", size=(truth_point_size*1.3), shape=17) +
  geom_point(data=df[df$type=="truth",], aes(x=x3, y=x4, color=indiv), size=truth_point_size, shape=17) +
  theme(legend.position = "none") +
  xlab("PCoA 3") +
  ylab("PCoA 4")
p <- grid.arrange(p1, p2, nrow=1)
if(setup == 1) {
  ggsave(paste0("output/same_baseline_diff_dynamics_PCA.png"), p, units="in", dpi=150, height=6, width=13)
} else {
  ggsave(paste0("output/diff_baseline_same_dynamics_PCA.png"), p, units="in", dpi=150, height=6, width=13)
}

plot_Sigma <- function(cov_mat, title) {
  temp <- gather_array(cov_mat, "value", "taxon1", "taxon2")
  p <- ggplot(temp, aes(taxon1, taxon2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low="darkblue", high="darkred", name="covariance") +
    theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(size = 10),
          text = element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
  return(p)
}

# compare Sigma and average Sigma_hat visually
plist_truth <- list()
plist_fitted <- list()
Kplus <- K
if(do_spike_in) {
  Kplus <- K+1
}
avg_Sigma <- matrix(NA, D-1, (D-1)*Kplus)
for(k in 1:Kplus) {
  # true Sigma
  truth_idx1 <- (k-1)*(D-1)*(n_samples+1) + 1
  truth_idx2 <- truth_idx1 + (D-1) - 1
  plist_truth[[k]] <- plot_Sigma(Sigmas[,truth_idx1:truth_idx2], title=paste0("truth (",k,")"))
  # average Sigma -- literally calculate the average over samples
  # avg_Sigma_hat <- matrix(0, D-1, D-1)
  # fit_idx1 <- truth_idx2 + 1
  # fit_idx2 <- fit_idx1 + (D-1) - 1
  # for(j in 1:n_samples) {
  #   Sigma_sample <- Sigmas[,fit_idx1:fit_idx2]
  #   avg_Sigma_hat <- avg_Sigma_hat + Sigma_sample
  #   fit_idx1 <- fit_idx2 + 1
  #   fit_idx2 <- fit_idx1 + (D-1) - 1
  # }
  # avg_Sigma_hat <- avg_Sigma_hat / n_samples
  # avg_Sigma[,((k-1)*(D-1)+1):(k*(D-1))] <- avg_Sigma_hat
  avg_Sigma_hat <- fits_MAP[[k]]$Sigma[,,1]
  avg_Sigma[,((k-1)*(D-1)+1):(k*(D-1))] <- avg_Sigma_hat
  plist_fitted[[k]] <- plot_Sigma(avg_Sigma_hat, title=paste0("fitted (",k,")"))
}

p <- do.call("grid.arrange", c(append(plist_truth, plist_fitted), ncol=(Kplus)))
if(setup == 1) {
  ggsave(paste0("output/same_baseline_diff_dynamics_Sigmas.png"), p, units="in", dpi=150, height=3.5, width=15)
} else {
  ggsave(paste0("output/diff_baseline_same_dynamics_Sigmas.png"), p, units="in", dpi=150, height=3.5, width=15)
}

# how tight are these posteriors? plot some estimates; eh... not really helpful

idx1 <- 3 # data set output k of K
idx2 <- 6 # data set output k' of K

# MAP estimates of Sigma (where the baseline is the same)
left <- (D-1)*(n_samples+1)*(idx1-1)+1
right <- left + (D-1) - 1
Sigma1 <- Sigmas[,left:right]
left <- (D-1)*(n_samples+1)*(idx2-1)+1
right <- left + (D-1) - 1
Sigma2 <- Sigmas[,left:right]

# posterior samples of Sigma
plist <- list()
plist[[length(plist)+1]] <- plot_Sigma(Sigma1, title="baseline")
plist[[length(plist)+1]] <- plot_Sigma(fits_MAP[[idx1]]$Sigma, title="MAP")
for(i in 1:10) {
  plist[[length(plist)+1]] <- plot_Sigma(fits[[idx1]]$Sigma[,,i], title="sample")
}
plist[[length(plist)+1]] <- plot_Sigma(Sigma1, title="baseline")
plist[[length(plist)+1]] <- plot_Sigma(fits_MAP[[idx2]]$Sigma, title="MAP")
for(i in 1:10) {
  plist[[length(plist)+1]] <- plot_Sigma(fits[[idx2]]$Sigma[,,i], title="sample")
}

p <- do.call("grid.arrange", c(plist, ncol=12))


# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
#                                                    DIAGNOSTICS
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------



if(FALSE & setup == 2) {

  idx1 <- 2  # pick a well estimated data set
  idx2 <- 10 # and a poorly estimated
  
  # -----------------------------------------------------------------------------------------------------------------
  #   DIAGNOSTICS: EVALUATE ETA
  # -----------------------------------------------------------------------------------------------------------------

  # MAP estimates of Eta
  par(mfrow=c(2,2))
  image(Etas[[idx1]])
  image(fits_MAP[[idx1]]$Eta[,,1])
  image(Etas[[idx2]])
  image(fits_MAP[[idx2]]$Eta[,,1])
  
  # posterior samples of Eta
  par(mfrow=c(4,3))
  image(Etas[[idx1]])
  for(i in 1:5) {
    image(fits[[idx1]]$Eta[,,i])
  }
  image(Etas[[idx2]])
  for(i in 1:5) {
    image(fits[[idx2]]$Eta[,,i])
  }
  
  # ordination of Eta
  eta_samples <- matrix(NA, 2*n_samples+2, (D-1)*N)
  eta_samples[1,] <- c(Etas[[idx1]])
  eta_labels <- c("baseline")
  for(i in 1:n_samples) {
    eta_samples[i+1,] <- c(fits[[idx1]]$Eta[,,i])
    eta_labels <- c(eta_labels, "fit1")
  }
  eta_samples[n_samples+2,] <- c(Etas[[idx2]])
  eta_labels <- c(eta_labels, "baseline")
  for(i in 1:n_samples) {
    eta_samples[i+n_samples+2,] <- c(fits[[idx2]]$Eta[,,i])
    eta_labels <- c(eta_labels, "fit2")
  }
  d_eta <- dist(eta_samples)
  embedding <- cmdscale(d_eta, k=2)
  df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                   label=as.factor(eta_labels))
  fitted_point_sz <- 2
  truth_point_size <- 4
  p1 <- ggplot() +
    geom_point(data=df[df$label!="baseline",], aes(x=x1, y=x2, color=label), size=fitted_point_sz) +
    geom_point(data=df[df$label=="baseline",], aes(x=x1, y=x2), size=truth_point_size, shape=17) +
    xlab("PCoA 1") +
    ylab("PCoA 2")
  p1
  
  # -----------------------------------------------------------------------------------------------------------------
  #   DIAGNOSTICS: EVALUATE LAMBDA
  # -----------------------------------------------------------------------------------------------------------------
  
  # MAP estimates of Lambda
  par(mfrow=c(2,2))
  image(Lambdas[[idx1]])
  image(fits_MAP[[idx1]]$Lambda[,,1])
  image(Lambdas[[idx2]])
  image(fits_MAP[[idx2]]$Lambda[,,1])
  
  # posterior samples of Eta
  par(mfrow=c(4,3))
  image(Lambdas[[idx1]])
  for(i in 1:5) {
    image(fits[[idx1]]$Lambda[,,i])
  }
  image(Lambdas[[idx2]])
  for(i in 1:5) {
    image(fits[[idx2]]$Lambda[,,i])
  }
  
  # ordination of Eta
  Lambda_samples <- matrix(NA, 2*n_samples+2, (D-1)*N)
  Lambda_samples[1,] <- c(Lambdas[[idx1]])
  Lambda_labels <- c("baseline")
  for(i in 1:n_samples) {
    Lambda_samples[i+1,] <- c(fits[[idx1]]$Lambda[,,i])
    Lambda_labels <- c(Lambda_labels, "fit1")
  }
  Lambda_samples[n_samples+2,] <- c(Lambdas[[idx2]])
  Lambda_labels <- c(Lambda_labels, "baseline")
  for(i in 1:n_samples) {
    Lambda_samples[i+n_samples+2,] <- c(fits[[idx2]]$Lambda[,,i])
    Lambda_labels <- c(Lambda_labels, "fit2")
  }
  d_Lambda <- dist(Lambda_samples)
  embedding <- cmdscale(d_Lambda, k=2)
  df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                   label=as.factor(Lambda_labels))
  fitted_point_sz <- 2
  truth_point_size <- 4
  p1 <- ggplot() +
    geom_point(data=df[df$label!="baseline",], aes(x=x1, y=x2, color=label), size=fitted_point_sz) +
    geom_point(data=df[df$label=="baseline",], aes(x=x1, y=x2), size=truth_point_size, shape=17) +
    xlab("PCoA 1") +
    ylab("PCoA 2")
  p1
  
  # -----------------------------------------------------------------------------------------------------------------
  #   DIAGNOSTICS: EVALUATE LAMBDA
  # -----------------------------------------------------------------------------------------------------------------
  
  # MAP estimates of Sigma (where the baseline is the same)
  left <- (D-1)*(n_samples+1)*(idx1-1)+1
  right <- left + (D-1) - 1
  Sigma1 <- Sigmas[,left:right]
  left <- (D-1)*(n_samples+1)*(idx2-1)+1
  right <- left + (D-1) - 1
  Sigma2 <- Sigmas[,left:right]
  
  par(mfrow=c(2,2))
  image(Sigma1)
  image(fits_MAP[[idx1]]$Sigma[,,1])
  image(Sigma2)
  image(fits_MAP[[idx2]]$Sigma[,,1])
  
  # posterior samples of Sigma
  par(mfrow=c(4,3))
  image(Sigma1)
  for(i in 1:5) {
    image(fits[[idx1]]$Sigma[,,i])
  }
  image(Sigma2)
  for(i in 1:5) {
    image(fits[[idx2]]$Sigma[,,i])
  }
  
  # ordination of Lambda
  Sigma_samples <- matrix(NA, 2*n_samples+2, (D-1)*(D-1))
  Sigma_samples[1,] <- c(Sigma1)
  Sigma_labels <- c("baseline")
  for(i in 1:n_samples) {
    Sigma_samples[i+1,] <- c(fits[[idx1]]$Sigma[,,i])
    Sigma_labels <- c(Sigma_labels, "fit1")
  }
  Sigma_samples[(n_samples+2),] <- c(Sigma2)
  Sigma_labels <- c(Sigma_labels, "baseline")
  for(i in 1:n_samples) {
    Sigma_samples[i+n_samples+2,] <- c(fits[[idx2]]$Sigma[,,i])
    Sigma_labels <- c(Sigma_labels, "fit2")
  }
  d_Sigma <- dist(Sigma_samples)
  embedding <- cmdscale(d_Sigma, k=2)
  df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                   label=as.factor(Sigma_labels))
  fitted_point_sz <- 2
  truth_point_size <- 4
  p1 <- ggplot() +
    geom_point(data=df[df$label!="baseline",], aes(x=x1, y=x2, color=label), size=fitted_point_sz) +
    geom_point(data=df[df$label=="baseline",], aes(x=x1, y=x2), size=truth_point_size, shape=17) +
    xlab("PCoA 1") +
    ylab("PCoA 2")
  p1
  
  # -----------------------------------------------------------------------------------------------------------------
  #   DIAGNOSTICS: WHICH TERMS IN SIGMA UPDATE CONTRIBUTE MORE TO VARIATION IN SIGMA?
  # -----------------------------------------------------------------------------------------------------------------
  
  # true
  
  temp <- Etas[[idx1]] - fits_MAP[[idx1]]$Lambda[,,1]
  upd1_1 <- temp%*%t(temp)

  temp <- fits_MAP[[idx1]]$Lambda[,,1] - Theta_fn(Xs[[idx1]])
  upd1_2 <- temp%*%solve(Gamma_fn(Xs[[idx1]]))%*%t(temp)
    
  temp <- Etas[[idx2]] - fits_MAP[[idx2]]$Lambda[,,1]
  upd2_1 <- temp%*%t(temp)
  
  temp <- fits_MAP[[idx2]]$Lambda[,,1] - Theta_fn(Xs[[idx2]])
  upd2_2 <- temp%*%solve(Gamma_fn(Xs[[idx2]]))%*%t(temp)

  par(mfrow=c(2,2))
  image(upd1_1)
  image(upd1_2)
  image(upd2_1)
  image(upd2_2)
  
  # -----------------------------------------------------------------------------------------------------------------
  #   END DIAGNOSTICS !!!
  # -----------------------------------------------------------------------------------------------------------------

}

# -----------------------------------------------------------------------------------------------------------------
#   OTHER EXPLORATORY STUFF: ARE THE DISTANCES BETWEEN SIGMA ESTIMATES RELATED TO DISTANCES IN BASELINE?
# -----------------------------------------------------------------------------------------------------------------

if(FALSE & setup == 2) {
  # (a) look at the distance matrices by eye: (L) baseline, (R) dynamics
  par(mfrow=c(1,2))
  d1 <- as.matrix(dist(t(Thetas)))
  image(d1)
  d2 <- matrix(NA, K, K)
  for(i in 1:K) {
    for(j in 1:K) {
      d2[i,j] <- Riemann_dist_pair(avg_Sigma[,((i-1)*(D-1)+1):(i*(D-1))], avg_Sigma[,((j-1)*(D-1)+1):(j*(D-1))])
    }
  }
  image(d2)

  # (b) does the baseline ordination resemble the dynamics ordination? (this is a long shot)
  # embed Theta
  d <- dist(t(Thetas)) # K x K
  embedding <- cmdscale(d, k=4)
  df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                   x3=embedding[,3], x4=embedding[,4],
                   indiv=as.factor(1:K))
  fitted_point_sz <- 2
  truth_point_size <- 4
  p1 <- ggplot(df) +
    geom_point(aes(x=x1, y=x2), color="black", size=(truth_point_size*1.3), shape=17) +
    geom_point(aes(x=x1, y=x2, color=indiv), size=truth_point_size, shape=17) +
    xlab("PCoA 1") +
    ylab("PCoA 2")
  p2 <- ggplot(df) +
    geom_point(aes(x=x3, y=x4), color="black", size=(truth_point_size*1.3), shape=17) +
    geom_point(aes(x=x3, y=x4, color=indiv), size=truth_point_size, shape=17) +
    theme(legend.position = "none") +
    xlab("PCoA 3") +
    ylab("PCoA 4")
  p <- grid.arrange(p1, p2, nrow=1)
  ggsave(paste0("output/diff_baseline_same_dynamics_PCA_Theta.png"), p, units="in", dpi=150, height=6, width=13)
}





