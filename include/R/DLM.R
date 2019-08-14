# ====================================================================================================================
# DYNAMIC LINEAR MODELS
# ====================================================================================================================

require("include.R")

# quick plot function for log ratios (x=time)
plot_lr <- function(ys, filename=NULL) {
  df <- gather_array(ys, "value", "time", "component")
  df$component <- as.factor(df$component)
  
  p <- ggplot(df, aes(x=time, y=value, group=component)) +
    geom_line(aes(color=component)) +
    theme_minimal() +
    theme(legend.position="none")
  if(!is.null(filename)) {
    ggsave(filename, scale=1.5, height=2, width=4)
  } else {
    show(p)
  }
}

# quick plot function for proportions (x=time)
plot_prop <- function(ys, filename=NULL) {
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
  if(!is.null(filename)) {
    ggsave(filename, scale=1.5, height=2, width=4)
  } else {
    show(p)
  }
}

# quick plot function for true thetas, filtered thetas, and smoothed thetas - sanity check
# data_obj is the output of five_taxa_simuation()
# fit_obj.f is the output of fit_filter()
# fit_obj.s is the output of fit_smoother()
plot_theta_fits <- function(data_obj, fit_obj.f=NULL, fit_obj.s=NULL, observation_vec=NULL, filename=NULL) {
  if(is.null(fit_obj.f) && is.null(fit_obj.s)) {
    return()
  }
  T <- nrow(data_obj$ys)
  T_actual <- T
  if(!is.null(observation_vec)) {
    T_actual <- max(observation_vec)
  }
  xs <- 1:T_actual
  D <- ncol(data_obj$ys)
  
  df.true <- gather_array(data_obj$ys, "logratio", "timepoint", "taxon")
  # replace timepoint with the actual day offset
  if(!is.null(observation_vec)) {
    for(i in 1:nrow(df.true)) {
      df.true$timepoint[i] <- observation_vec[df.true$timepoint[i]]
    }
  }
  df.true <- cbind(df.true, which="true")
  
  filtered.observations <- matrix(0, T_actual, D)
  smoothed.observations <- matrix(0, T_actual, D)
  for(t in 1:T_actual) {
    if(!is.null(fit_obj.f)) {
      filtered.observations[t,] <- data_obj$F%*%fit_obj.f$Thetas.t[,,t]
    }
    if(!is.null(fit_obj.s)) {
      smoothed.observations[t,] <- data_obj$F%*%fit_obj.s$Thetas.t[,,t]
    }
  }
  df.all <- df.true
  if(!is.null(fit_obj.f)) {
    df.filtered <- gather_array(filtered.observations, "logratio", "timepoint", "taxon")
    df.filtered <- cbind(df.filtered, which="filtered")
    df.all <- rbind(df.all, df.filtered)
  }
  if(!is.null(fit_obj.s)) {
    df.smoothed <- gather_array(smoothed.observations, "logratio", "timepoint", "taxon")
    df.smoothed <- cbind(df.smoothed, which="smoothed")
    df.all <- rbind(df.all, df.smoothed)
  }
  
  p <- ggplot(data=df.all, aes(x=timepoint, y=logratio, color=which, group=which)) + 
    geom_line() + 
    facet_wrap(~taxon) +
    theme_minimal()
  if(!is.null(filename)) {
    ggsave(filename, plot=p, scale=1.5, width=8, height=5)
  } else {
    show(p)
  }
}

# unfinished; busted
# fit_and_plot_intervals <- function(data_obj, observation_vec=observation_vec) {
#   S <- 10
#   T <- nrow(data_obj$ys)
#   if(!is.null(observation_vec)) {
#     T <- max(observation_vec)
#   }
#   draws <- array(dim=c(T, ncol(data_obj$ys), S)) # time (rows) x taxa (columns) x samples
#   for(s in 1:S) {
#     cat("Generating sample",s,"\n")
#     fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
#     fit.s <- fit_smoother(data_obj, fit.f)
#     draws[,,s] <- t(fit.s$Thetas.t[1,,])
#   }
#   labeled_ys <- data_obj$ys
#   if(!is.null(observation_vec)) {
#     rownames(labeled_ys) <- observation_vec
#   }
#   samples <- gather_array(data_obj$ys, "coord", "timepoint", "taxon")
#   # crazy inefficient; fix this
#   samples <- cbind(samples, lower=NA)
#   samples <- cbind(samples, upper=NA)
#   for(i in 1:nrow(samples)) {
#     samples[i,"lower"] <- min(draws[samples[i,"timepoint"],samples[i,"taxon"],])
#     samples[i,"upper"] <- max(draws[samples[i,"timepoint"],samples[i,"taxon"],])
#   }
#   taxon <- 2
#   p <- ggplot(samples[samples$taxon==taxon,], aes(timepoint)) +
#     geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey80") +
#     geom_line(aes(y=coord))
#   p
# }

plot_mean_Sigma <- function(fit_obj, filename=NULL) {
  mean.Sigma <- fit_obj$Xi/(fit_obj$upsilon - nrow(fit_obj$Xi) - 1)
  if(!is.null(filename)) {
    png(filename)
    image(mean.Sigma)
    dev.off()
  } else {
    image(mean.Sigma)
  }
}


# generate Fourier form rotation matrix for a given omega (a function of period)
build_G <- function(period, harmonics=1) {
  if(harmonics == 1) {
    omega <- 2*pi/period
    return(matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2))
  } else {
    G <- matrix(0, 2*harmonics, 2*harmonics)
    for(h in 1:harmonics) {
      omega <- 2*pi*h/period
      common_offset <- (h-1)*2
      G[(common_offset+1),(common_offset+1)] <- cos(omega)
      G[(common_offset+1),(common_offset+2)] <- sin(omega)
      G[(common_offset+2),(common_offset+1)] <- -sin(omega)
      G[(common_offset+2),(common_offset+2)] <- cos(omega)
    }
    return(G)
  }
}

# 2-taxa simulation just for illustrative purposes
two_taxa_simulation <- function(T=30) {
  state_noise_scale <- 0.1
  gamma.t <- 1
  observational_noise_scale <- state_noise_scale*gamma.t
  
  period <- 10

  F <- matrix(c(1, 0), 1, 2)
  W.t <- diag(2)*state_noise_scale
  G <- build_G(period)
  upsilon <- 100
  
  # we'll do all combinations of noise and generate plots to see how the resultant log ratios
  # and proportions look
  noise_flags <- matrix(FALSE, 4, 2)
  noise_flags[2,1] <- TRUE # use state transition noise
  noise_flags[3,2] <- TRUE # use observational noise
  noise_flags[4,] <- TRUE # use both types of noise
  
  D <- 2
  
  for(i in 1:dim(noise_flags)[1]) {
    Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
    Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
    
    theta.t <- rmatrixnormal(1, matrix(1, 2, D), W.t, Sigma)[,,1]

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

five_taxa_simulation <- function(T=40, indep_taxa=TRUE, uniform_start=TRUE, noise_scale=0.05, save_images=TRUE) {
  state_noise_scale <- noise_scale
  gamma <- 1
  observational_noise_scale <- state_noise_scale*gamma
  
  period <- 10

  F <- matrix(c(1, 0), 1, 2)
  W <- diag(2)*state_noise_scale
  G <- build_G(period)
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
  
  # very similar initial states; everybody is pretty much in phase here
  if(uniform_start) {
    M.0 <- matrix(1, 2, D)
    tag <- paste0(tag, "unifinit")
  } else {
    M.0 <- matrix(rnorm(10, 1, 1), 2, D)
    tag <- paste0(tag, "nonunifinit")
  }
  C.0 <- W
  
  theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1]
  
  ys <- matrix(0, T, D)
  for(t in 1:T) {
    # state equation
    theta.t <- G%*%theta.t + rmatrixnormal(1, matrix(0, 2, D), W, Sigma)[,,1]
    # observation equation
    ys[t,] <- F%*%theta.t + rmvnorm(1, rep(0, D), gamma*Sigma)
  }
  
  if(save_images) {
    plot_lr(ys, filename=paste0("plots/DLMsim_logratios_",tag,".png"))
    plot_prop(ys, filename=paste0("plots/DLMsim_proportions_",tag,".png"))
  }
  
  return(list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=Sigma, M.0=M.0, C.0=C.0))
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

build_A <- function(T, G, C.0, W, gamma, save_images=T) {
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
            addend <- G_left%*%W%*%G_right
            first_sum <- first_sum + addend
          }
        }
        # second sum
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, t, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C.0%*%G_right
        A[t,t] <- gamma + F%*%(W + first_sum + second_sum)%*%t(F)
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
          first_sum <- first_sum + G_left%*%W%*%G_right
        }
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, tk, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C.0%*%G_right
        G_left <- stack_G(G, t, tk+1)
        A[i,j] <- F%*%(G_left%*%W + first_sum + second_sum)%*%t(F)
      }
    }
  }
  
  if(save_images) {
    png("plots/DLMsim_Amat.png")
    image(A)
    dev.off()
  }
  
  if(min(eigen(A)$values) < 0) {
    cat("Matrix A has negative eigenvalue(s)!\n")
  }
}

# Kalman filter
# data_obj is a list containing ys, F, W, G, upsilon, Xi, gamma, Sigma, M.0, C.0
# observation_vec (if present) indicates the spacing of observations, e.g. c(1, 3, 4, 7)
#   indicates observations 2, 5, & 6 are missing and should be imputed in the usual way
fit_filter <- function(data_obj, censor_vec=NULL, observation_vec=NULL, discount=NULL) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  if(!is.null(observation_vec)) {
    T <- max(observation_vec)
  } else {
    T <- nrow(data_obj$ys)
  }
  if(is.null(censor_vec)) {
    censor_vec <- rep(0, T)
  }
  upsilon.t <- data_obj$upsilon
  Xi.t <- data_obj$Xi
  M.t <- data_obj$M.0
  C.t <- data_obj$C.0
  Thetas.t <- array(0, dim=c(Theta.dim, D, T)) # sample at each t
  Cs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  Rs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  W.t <- NULL # last observed W.t (when using discount)
  for(t in 1:T) {
    if(censor_vec[t] == 1 || (!is.null(observation_vec) && !(t %in% observation_vec))) {
      # case 1: no observations at this time point for any individuals
      # impute all by passing prior on as posterior
      P.t <- data_obj$G%*%C.t%*%t(data_obj$G)
      if(is.null(discount)) {
        R.t <- P.t + data_obj$W
      } else {
        # note: West & Harrison (p. 352) suggest not applying discount for missing observations
        if(is.null(W.t)) {
          W.t <- ((1-discount)/discount)*P.t
        }
        R.t <- P.t + W.t
      }
      R.t <- round(R.t, 10)
      Rs.t[,,t] <- R.t
      A.t <- data_obj$G%*%M.t
      # carry forward without update
      M.t <- A.t
      Ms.t[,,t] <- M.t
      C.t <- R.t
      Cs.t[,,t] <- C.t
      # no change to Sigma.t parameters
      Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
      Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
    } else {
      # case 2: apply update where possible, note: F.t.T is F for us here and G.t is G
      # one-step ahead forecast at t
      P.t <- data_obj$G%*%C.t%*%t(data_obj$G)
      if(is.null(discount)) {
        R.t <- P.t + data_obj$W
      } else {
        W.t <- ((1-discount)/discount)*P.t
        R.t <- P.t + W.t
      }
      R.t <- round(R.t, 10)
      Rs.t[,,t] <- R.t
      A.t <- data_obj$G%*%M.t
      f.t.T <- data_obj$F%*%A.t
      q.t <- data_obj$gamma + (data_obj$F%*%R.t%*%t(data_obj$F))[1,1]
      if(is.null(observation_vec)) {
        # we're assuming 1 individual if observation_mat is missing
        e.t.T <- data_obj$ys[t,] - f.t.T
      } else {
        # if more than one sample on the same day (mislabeled?) take the first for now
        e.t.T <- data_obj$ys[as(which(observation_vec == t)[1], "numeric"),] - f.t.T
      }
      S.t <- R.t%*%t(data_obj$F)/q.t
      M.t <- A.t + S.t%*%e.t.T
      Ms.t[,,t] <- M.t
      C.t <- R.t - q.t*S.t%*%t(S.t)
      C.t <- round(C.t, 10)
      Cs.t[,,t] <- C.t
      upsilon.t <- upsilon.t + 1
      Xi.t <- Xi.t + t(e.t.T)%*%e.t.T/q.t
      Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
      Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
    }
  }
  return(list(Thetas.t=Thetas.t, upsilon=upsilon.t, Xi=Xi.t, Ms.t=Ms.t, Cs.t=Cs.t, Rs.t=Rs.t))
}

# simulation smoother
# fit_obj is the output of fit_filter()
fit_smoother <- function(data_obj, fit_obj) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  T <- dim(fit_obj$Thetas.t)[3]
  Sigma.t <- rinvwishart(1, fit_obj$upsilon, fit_obj$Xi)[,,1]
  Thetas.t.smoothed <- array(0, dim=c(Theta.dim, D, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  etas.t <- matrix(0, D, T)
  rmatnorm_mean <- fit_obj$Ms.t[,,T]
  dim(rmatnorm_mean) <- c(Theta.dim, D) # fix if the state has one dimension only
                                        # can't do drop=FALSE here
  Thetas.t.smoothed[,,T] <- rmatrixnormal(1, rmatnorm_mean, fit_obj$Cs.t[,,T], Sigma.t)[,,1]
  etas.t[,T] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,T], 1, Sigma.t)[,,1]
  Ms.t[,,T] <- rmatnorm_mean
  for(t in (T-1):1) {
    Z.t <- fit_obj$Cs.t[,,t]%*%t(data_obj$G)%*%solve(fit_obj$Rs.t[,,(t+1)])
    M.t.star <- fit_obj$Ms.t[,,t] + Z.t%*%(Thetas.t.smoothed[,,(t+1)] - data_obj$G%*%fit_obj$Ms.t[,,t])
    Ms.t[,,t] <- M.t.star
    C.t.star <- round(fit_obj$Cs.t[,,t] - Z.t%*%fit_obj$Rs.t[,,(t+1)]%*%t(Z.t), 10)
    Thetas.t.smoothed[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.t)[,,1]
    # draw an eta, this is how we'll estimate modeled covariance
    etas.t[,t] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,t], 1, Sigma.t)[,,1]
  }
  return(list(Thetas.t=Thetas.t.smoothed, etas.t=etas.t, Ms.t=Ms.t))
}

pull_indiv_data <- function(sname, asv_data, pseudocount=0.5, date_begin="1900-01-01", date_end="2100-01-01",
                               subset_dim=0, lr_transform=TRUE) {
  pruned <- prune_samples(sample_data(non_reps)$sname==sname, non_reps)
  pruned <- prune_samples((sample_data(pruned)$collection_date > date_begin) &
                            (sample_data(pruned)$collection_date < date_end), pruned)
  cat(paste0("Real data set (",sname,") has ",nsamples(pruned)," samples and ",ntaxa(pruned)," taxa\n"))
  counts <- otu_table(pruned)@.Data # samples (rows) x taxa (columns)
  if(lr_transform) {
    ys <- driver::clr(counts + pseudocount)
  } else {
    ys <- counts
  }
  # the baseline date common to individuals
  min_date <- min(sample_data(pruned)$collection_date)
  dates_observed <- sample_data(pruned)$collection_date
  observation_vec <- sapply(dates_observed, function(x) { round(difftime(x, min_date, units="days"))+1 } )
  # subset taxa
  if(subset_dim > 0) {
    ys <- ys[,1:subset_dim]
  }
  return(list(ys=ys, observation_vec=observation_vec, sname=sname))
}
