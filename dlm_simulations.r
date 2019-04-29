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

state_noise_scale <- 0.1
gamma.t <- 1
observational_noise_scale <- state_noise_scale*gamma.t

period <- 10
omega <- 2*pi*1/period

F <- matrix(c(1, 0), 1, 2)
G <- build_G(omega)
W.t <- diag(2)*state_noise_scale
upsilon <- 100

# two taxa - very simple example
if(FALSE) {
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

# five taxa, longer simulation
if(TRUE) {
  state_noise <- TRUE
  observational_noise <- TRUE
  
  state_noise_scale <- 0.05
  gamma.t <- 1
  observational_noise_scale <- state_noise_scale*gamma.t
  
  period <- 10
  omega <- 2*pi*1/period
  
  F <- matrix(c(1, 0), 1, 2)
  G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)
  W.t <- diag(2)*state_noise_scale
  upsilon <- 100
  
  D <- 5
  
  Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  Sigma.true <- Sigma
  
  # very similar initial states; everybody is pretty much in phase here
  M.0 <- matrix(1, 2, D)
  C.0 <- W.t
  
  # different initial states (random)
  #C.0 <- W.t*100
  
  theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1]

  T <- 30
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
  
  plot_lr(ys, filename="big_1.png")
  
  plot_prop(ys, filename="big_2.png")
}

stack_G_old <- function(G, it_begin, it_end, descending=TRUE, transpose=FALSE) {
  obj <- G
  if(transpose) {
    obj <- t(G)
  }
  if(it_begin != it_end) {
    if(descending) {
      if(it_begin > it_end) {
        ret_val <- obj
        for(i in 1:(it_begin-it_end)) {
          ret_val <- ret_val%*%obj
        }
      } else {
        # invalid case
        ret_val <- matrix(0, 2, 2)
      }
    } else {
      if(it_begin < it_end) {
        ret_val <- obj
        for(i in 1:(it_end-it_begin)) {
          ret_val <- ret_val%*%obj
        }
      } else {
        # invalid case
        ret_val <- matrix(0, 2, 2)
      }
    }
  } else {
    ret_val <- obj
  }
  #ret_val <- Re(e_vec%*%e_val%*%t(e_vec))
  return(ret_val)
}

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
  ret_val <- Re(e_vec%*%e_val%*%solve(e_vec))
  return(ret_val)
}

k_per <- function(x, xp) {
  sigma.sq <- 1
  ell.sq <- 1
  period <- 10
  euclid_dist <- dist(rbind(x, xp))
  return(sigma.sq * exp(-(2*sin(pi * euclid_dist / period)**2)/ell.sq))
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
          G_left <- stack_G(G, t, 2)
          G_right <- stack_G(G, ell, t, descending=FALSE, transpose=TRUE)
          first_sum <- first_sum + G_left%*%W.t%*%G_right
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

df <- data.frame(x=seq(1,T), y=A[1,], time="t1")
df <- rbind(df, data.frame(x=seq(1,T), y=A[3,], time="t3"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[5,], time="t5"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[7,], time="t7"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[10,], time="t10"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[12,], time="t12"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[15,], time="t15"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[17,], time="t17"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[20,], time="t20"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[25,], time="t25"))
df <- rbind(df, data.frame(x=seq(1,T), y=A[30,], time="t30")) # what's going on with this time point
p <- ggplot(df, aes(x=x, y=y, color=time)) +
  geom_line() + 
  theme_minimal()
p
ggsave(paste("covariance_entries_DLM.png", sep=""), width=6, height=4, units="in", scale=1.5)
eigen(A)$values

# check for diagnonal dominance - HELL NO!
for(i in 1:T) {
  diag <- abs(A[i,i])
  off_diag <- 0
  for(j in 1:T) {
    if(j != i) {
      off_diag <- off_diag + abs(A[i,j])
    }
  }
  cat("Diag:",diag,", off-diag:",off_diag,"\n")
}

dd_fail <- 0
# check diagonal dominance
for(i in 1:T) {
  off_diag <- sum(A[i,]) - A[i,i]
  if(off_diag > A[i,i]) {
    dd_fail <- dd_fail + 1
  }
}
cat("Diagonal fails to dominate for:",dd_fail,"\n")

# periodic kernel
A2 <- matrix(0, T, T)
for(i in 1:T) {
  for(j in 1:T) {
    A2[i,j] <- k_per(ys[i,], ys[j,])
  }
}
image(A2)
eigen(A2)$values
is.positive.semidefinite(A2)

# Kalman filter for Sigma
# empirical covariance
Xi.t <- t(ys)%*%ys/T
upsilon.t <- D + 2

# THIS APPEARS TO BE WRONG, CHECK IT!
M.tm1 <- M.0
C.tm1 <- C.0
for(t in 1:T) {
  # note: F.t.T is F for us here and G.t is G
  # prior at t
  A.t <- G%*%M.tm1
  R.t <- G%*%C.tm1%*%t(G) + W.t
  Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
  theta.t <- rmatrixnormal(1, A.t, R.t, Sigma)[,,1]
  # one-step ahead forecast at t
  f.t.T <- F%*%A.t
  q.t <- gamma.t + (F%*%R.t%*%t(F))[1,1]
  # posterior at t
  e.t.T <- ys[t,] - f.t.T
  S.t <- R.t%*%t(F)/q.t
  M.t <- A.t + S.t%*%e.t.T
  C.t <- R.t - q.t*S.t%*%t(S.t)
  upsilon.t <- upsilon.t + 1
  Xi.t <- Xi.t + t(e.t.T)%*%e.t.T/q.t
}
Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
image(Sigma.true)
image(Sigma.t)


















