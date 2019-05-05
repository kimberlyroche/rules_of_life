library(LaplacesDemon)
library(mvtnorm)
library(matrixsampling)
library(driver)
library(ggplot2)
library(RColorBrewer)

build_G <- function(angle) {
  return(matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2))
}

plot_lr <- function(ys, filename) {
  df <- gather_array(ys, "value", "time", "component")
  df$component <- as.factor(df$component)
  
  p <- ggplot(df, aes(x=time, y=value, group=component)) +
    geom_line(aes(color=component)) +
    theme_minimal() +
    theme(legend.position="none")
  ggsave(filename, scale=1.5, height=2, width=4)
}

plot_prop <- function(ys, filename) {
  # visualize as proportions
  proportions <- clrInv(ys)
  df <- gather_array(proportions, "value", "time", "component")
  df$component <- factor(df$component)
  
  categories <- unique(df$component)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df$component)))
  p <- ggplot(df, aes(x=time, y=value, fill=component)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none")
  ggsave(filename, scale=1.5, height=2, width=4)
}

# two taxa - very simple example
if(FALSE) {
  state_noise_scale <- 0.1
  gamma.t <- 1
  observational_noise_scale <- state_noise_scale*gamma.t
  
  period <- 10
  omega <- 2*pi*1/period
  
  F <- matrix(c(1, 0), 1, 2)
  W.t <- diag(2)*state_noise_scale
  G <- build_G(omega)
  upsilon <- 100
  
  noise_flags <- matrix(FALSE, 4, 2)
  noise_flags[2,1] <- TRUE # use state transition noise
  noise_flags[3,2] <- TRUE # use observational noise
  noise_flags[4,] <- TRUE # use both types of noise
  
  D <- 2
  
  for(i in 1:dim(noise_flags)[1]) {
    Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
    Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
    
    theta.t <- rmatrixnormal(1, matrix(1, 2, D), W.t, Sigma)[,,1]
    
    T <- 30
    ys <- matrix(0, T, D)
    for(t in 1:T) {
      # state equation
      theta.t <- G%*%theta.t
      if(noise_flags[i,1]) {
        theta.t = theta.t + rmatrixnormal(1, matrix(0, 2, D), W.t, Sigma)[,,1]
      }
      # observation equation
      ys[t,] <- F%*%theta.t
      if(noise_flags[i,2]) {
        ys[t,] <- ys[t,] + rmvnorm(1, rep(0, D), gamma.t*Sigma)
      }
    }
    
    plot_lr(ys, filename=paste(i,"_1.png",sep=""))
    
    plot_prop(ys, filename=paste(i,"_2.png",sep=""))
  }
}

# five taxa, long simulation
if(TRUE) {
  indep_taxa <- TRUE
  uniform_start <- FALSE
  
  state_noise <- TRUE
  observational_noise <- TRUE
  
  state_noise_scale <- 0.05
  gamma.t <- 1
  observational_noise_scale <- state_noise_scale*gamma.t
  
  period <- 10
  omega <- 2*pi*1/period
  
  F <- matrix(c(1, 0), 1, 2)
  W.t <- diag(2)*state_noise_scale
  G <- build_G(omega)
  upsilon <- 100
  
  D <- 5
  
  if(indep_taxa) {
    # fully independent taxa
    Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
    tag <- "indeptax_"
  } else {
    # correlated/anti-correlated taxa
    Xi <- diag(D)
    Xi[1,2] <- 0.5
    Xi[2,1] <- Xi[1,2]
    Xi[3,4] <- -0.5
    Xi[4,3] <- Xi[3,4]
    Xi <- Xi*observational_noise_scale*(upsilon-D-1)
    tag <- "nonindeptax_"
  }
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  Sigma.true <- Sigma
  
  # very similar initial states; everybody is pretty much in phase here
  if(uniform_start) {
    M.0 <- matrix(1, 2, D)
    tag <- paste0(tag, "unifinit")
  } else {
    M.0 <- matrix(rnorm(10, 1, 1), 2, D)
    tag <- paste0(tag, "nonunifinit")
  }
  C.0 <- W.t
  
  theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1]

  T <- 40
  ys <- matrix(0, T, D)
  for(t in 1:T) {
    # state equation
    theta.t <- G%*%theta.t
    if(state_noise) {
      theta.t = theta.t + rmatrixnormal(1, matrix(0, 2, D), W.t, Sigma)[,,1]
    }
    # observation equation
    ys[t,] <- F%*%theta.t
    if(observational_noise) {
      ys[t,] <- ys[t,] + rmvnorm(1, rep(0, D), gamma.t*Sigma)
    }
  }
  
  plot_lr(ys, filename=paste0("plots/DLMsim_logratios_",tag,".png"))
  
  plot_prop(ys, filename=paste0("plots/DLMsim_proportions_",tag,".png"))
}

# powers G through eigenvalue decomposition
stack_G <- function(G, it_begin, it_end, descending=TRUE, transpose=FALSE) {
  obj <- G
  if(transpose) {
    obj <- t(G)
  }
  obj <- eigen(obj)
  e_vec <- obj$vectors
  e_val <- diag(2)
  diag(e_val) <- obj$values
  if(it_begin != it_end) {
    if(descending) {
      if(it_begin > it_end) {
        power_it <- length(1:(it_begin-it_end+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    } else {
      if(it_begin < it_end) {
        power_it <- length(1:(it_end-it_begin+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    }
  }
  # explicitly only returning the real part of A
  # some tiny complex eigenvalues can be produced -- cool to truncate in this way?
  ret_val <- Re(e_vec%*%e_val%*%solve(e_vec))
  return(ret_val)
}

# calculate A (covariance matrix over states) for this simulation
# this is the expression exactly as in the manuscript, calculated from Cov(eta_t, eta_{t-k})
A <- matrix(0, T, T)
for(i in 1:T) {
  for(j in 1:T) {
    if(i == j) {
      # diagonal
      t <- j
      first_sum <- matrix(0, 2, 2)
      if(t >= 2) {
        for(ell in t:2) {
          G_left <- stack_G(G, t, ell)
          G_right <- stack_G(G, ell, t, descending=FALSE, transpose=TRUE)
          addend <- G_left%*%W.t%*%G_right
          first_sum <- first_sum + addend
        }
      }
      # second sum
      G_left <- stack_G(G, t, 1)
      G_right <- stack_G(G, 1, t, descending=FALSE, transpose=TRUE)
      second_sum <- G_left%*%C.0%*%G_right
      A[t,t] <- gamma.t + F%*%(W.t + first_sum + second_sum)%*%t(F)
    } else {
      tk <- i
      t <- j
      if(j < i) {
        tk <- j
        t <- i
      }
      # off-diagonal
      first_sum <- matrix(0, 2, 2)
      for(ell in tk:2) {
        G_left <- stack_G(G, t, ell)
        G_right <- stack_G(G, ell, tk, descending=FALSE, transpose=TRUE)
        first_sum <- first_sum + G_left%*%W.t%*%G_right
      }
      G_left <- stack_G(G, t, 1)
      G_right <- stack_G(G, 1, tk, descending=FALSE, transpose=TRUE)
      second_sum <- G_left%*%C.0%*%G_right
      G_left <- stack_G(G, t, tk+1)
      A[i,j] <- F%*%(G_left%*%W.t + first_sum + second_sum)%*%t(F)
    }
  }
}

png("plots/DLMsim_Amat.png")
image(A)
dev.off()

eigen(A)$values

use_perfect_priors <- FALSE
if(use_perfect_priors){
  upsilon.t <- upsilon
  Xi.t <- Xi
  M.t <- M.0
  C.t <- W.t
} else {
  # need a good way of setting these priors
  upsilon.t <- D+2
  empirical_cov <- t(ys)%*%ys/T
  Xi.t <- empirical_cov*(upsilon.t-D-1) # center loosely at empirical covariance
  M.t <- M.0
  C.t <- diag(2)*mean(diag(empirical_cov))
}
Thetas.t <- array(0, dim=c(2, D, T)) # sample at each t
Cs.t <- array(0, dim=c(2, 2, T))
Ms.t <- array(0, dim=c(2, D, T))
Rs.t <- array(0, dim=c(2, 2, T))
for(t in 1:T) {
  # note: F.t.T is F for us here and G.t is G
  # prior at t
  A.t <- G%*%M.t
  R.t <- G%*%C.t%*%t(G) + W.t
  Rs.t[,,t] <- R.t
  # one-step ahead forecast at t
  f.t.T <- F%*%A.t
  q.t <- gamma.t + (F%*%R.t%*%t(F))[1,1]
  # posterior at t
  e.t.T <- ys[t,] - f.t.T
  S.t <- R.t%*%t(F)/q.t
  M.t <- A.t + S.t%*%e.t.T
  Ms.t[,,t] <- M.t
  C.t <- R.t - q.t*S.t%*%t(S.t)
  Cs.t[,,t] <- C.t
  upsilon.t <- upsilon.t + 1
  Xi.t <- Xi.t + t(e.t.T)%*%e.t.T/q.t
  Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
  Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
}

cat("Trace true Sigma:",sum(diag(Sigma.true)),"\n")
png(paste0("plots/DLMsim_true_Sigma_",tag,".png"), width=500, height=500)
image(Sigma.true)
dev.off()

mean.Sigma.t <- Xi.t/(upsilon.t - D - 1)
cat("Trace inferred Sigma:",sum(diag(mean.Sigma.t)),"\n")
png(paste0("plots/DLMsim_mean_Sigma_",tag,".png"), width=500, height=500)
image(mean.Sigma.t)
dev.off()

# simulation smoother
Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
Thetas.t.smoothed <- array(0, dim=c(2, D, T))
Thetas.t.smoothed[,,T] <- rmatrixnormal(1, Ms.t[,,T], Cs.t[,,T], Sigma)[,,1]
for(t in (T-1):1) {
  Z.t <- Cs.t[,,t]%*%t(G)%*%solve(Rs.t[,,(t+1)])
  M.t.star <- Ms.t[,,t] + Z.t%*%(Thetas.t.smoothed[,,(t+1)] - G%*%Ms.t[,,t])
  C.t.star <- round(Cs.t[,,t] - Z.t%*%Rs.t[,,(t+1)]%*%t(Z.t), 10)
  Thetas.t.smoothed[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.t)[,,1]
}

# evaluate the fitted thetas as a sanity check
df.true <- gather_array(ys, "logratio", "timepoint", "taxon")
df.true <- cbind(df.true, which="true")

filtered.observations <- matrix(0, T, D)
smoothed.observations <- matrix(0, T, D)
for(t in 1:T) {
  filtered.observations[t,] <- F%*%Thetas.t[,,t]
  smoothed.observations[t,] <- F%*%Thetas.t.smoothed[,,t]
}
df.filtered <- gather_array(filtered.observations, "logratio", "timepoint", "taxon")
df.filtered <- cbind(df.filtered, which="filtered")
df.smoothed <- gather_array(smoothed.observations, "logratio", "timepoint", "taxon")
df.smoothed <- cbind(df.smoothed, which="smoothed")

df.all <- rbind(df.true, df.filtered, df.smoothed)

p <- ggplot(data=df.all, aes(x=timepoint, y=logratio, color=which, group=which)) + 
  geom_line() +
  facet_wrap(~taxon) +
  theme_minimal()
p
ggsave(paste0("plots/DLMsim_theta_fits_",tag,".png"), plot=p, scale=1.5, width=8, height=5)















