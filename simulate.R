source("sim_densities.R")

# SIMULATION FUNCTIONS

simulate_univariate_normal <- function(it) {
  N <- it
  data <- rnorm(N, 0, 2)

  avg_est <- 0
  for(i in 1:10) {
    res <- optim(par=c(runif(1)), fn=logd_uni_normal, data=data, N=N, method="L-BFGS-B")
    avg_est <- avg_est + res$par
  }
  cat("Actual: 2\n")
  cat("Average estimate:",(avg_est/10),"\n")
}

simulate_multivariate_normal <- function(it) {
  N <- it
  D <- 5
  Sigma_true <- 2*diag(D)
  data <- t(rmvnorm(N, mean=rep(0,D), sigma=Sigma_true))

  avg_est <- 0
  for(i in 1:10) {
    res <- optim(par=c(runif(1)), fn=logd_mult_normal, data=data, N=N, D=D, method="L-BFGS-B")
    avg_est <- avg_est + res$par
  }
  cat("Actual: 2\n")
  cat("Average estimate:",(avg_est/10),"\n")
}

simulate_matrix_normal <- function(it) {
  P <- 20
  N <- 30
  zero_mean <- matrix(0, P, N)

  sigma_true.1 <- 0.25
  component.1 <- diag(N)
  sigma_true.2 <- 2
  component.2 <- matrix(0, nrow=N, ncol=N)
  component.2[1:15,1:15] <- 0.8
  component.2[16:30,16:30] <- 0.8
  component.2 <- component.2 + diag(N)*0.1

  sample_data <- array(0, dim=c(P, N, it))
  sample_covar <- round((sigma_true.1*component.1 + sigma_true.2*component.2), digits=10)
  for(i in 1:it) {
    sample_data[,,i] <- rmatrixnorm(zero_mean, diag(P), sample_covar)
  }

  avg_est.1 <- 0
  avg_est.2 <- 0
  for(i in 1:1) {
    cat("Opt it:",i,"\n")
    res <- optim(par=c(runif(1), runif(1)), fn=logd_mat_norm, vc=list(component.1, component.2), data=sample_data, P=P, method="L-BFGS-B")
    avg_est.1 <- avg_est.1 + res$par[1]
    avg_est.2 <- avg_est.2 + res$par[2]
  }
  cat("Actual: 0.25\n")
  cat("Average estimate:",exp(avg_est.1/1),"\n")
  cat("Actual: 2\n")
  cat("Average estimate:",exp(avg_est.2/1),"\n")
}

simulate_matT_one <- function(it) {
  N <- 20
  P <- 15
  upsilon <- P + 100
  zero_mean <- matrix(0, P, N)
  sigma_true.1 <- 2
  component.1 <- diag(N)
  #component.1[1:10,1:10] <- component.1[1:10,1:10] + 0.8
  #component.1[11:20,11:20] <- component.1[11:20,11:20] + 0.8
  sample_data <- array(0, dim=c(P, N, it))
  sample_covar.2 <- round((sigma_true.1*component.1), digits=10)
  for(i in 1:it) {
    sample_covar.1 <- round(riwish(upsilon, diag(P)), digits=10)
    sample_data[,,i] <- rmatrixnorm(zero_mean, sample_covar.1, sample_covar.2)
  }

  avg_est <- 0
  for(i in 1:1) {
    cat("Opt it:",i,"\n")
    res <- optim(par=c(runif(1)), fn=logd_mat_t_one, vc=list(component.1), data=sample_data, N=N, P=P, upsilon=upsilon, method="L-BFGS-B")
    avg_est <- avg_est + res$par
  }
  cat("Actual: 2\n")
  cat("Average estimate:",exp(avg_est/1),"\n")
}

simulate_matT_two <- function(it) {
  #N <- 20
  #P <- 15
  N <- 300
  P <- 20
  upsilon <- P + 2
  zero_mean <- matrix(0, P, N)
  # component 1 - AR1 date distance
  sigma_true.1 <- 0.2
  component.1 <- diag(N)
  for(i in 1:(N-1)) {
    component.1[i,i+1] <- 0.33
    component.1[i+1,i] <- 0.33
  }
  # component 1 - seasonal distance
  sigma_true.2 <- 0.8
  component.2 <- matrix(-0.2, nrow=N, ncol=N)
  component.2[1:150,1:150] <- 0.2
  component.2[151:300,151:300] <- 0.2
  for(i in 1:N) {
    component.2[i,i] <- 1
  }
  sample_data <- array(0, dim=c(P, N, it))
  sample_covar.2 <- round((sigma_true.1*component.1 + sigma_true.2*component.2), digits=5)
  for(i in 1:it) {
    sample_covar.1 <- round(riwish(upsilon, diag(P)), digits=5)
    sample_data[,,i] <- rmatrixnorm(zero_mean, sample_covar.1, sample_covar.2)
  }

  avg_est.1 <- 0
  avg_est.2 <- 0
  for(i in 1:1) {
    res <- optim(par=c(runif(1), runif(1)),
                 fn=logd_mat_t_two,
                 vc=list(component.1, component.2),
                 data=sample_data, N=N, P=P, upsilon=upsilon)
    avg_est.1 <- avg_est.1 + res$par[1]
    avg_est.2 <- avg_est.2 + res$par[2]
  }
  cat("Actual:",sigma_true.1,", Average:",exp(avg_est.1),"\n")
  cat("Actual:",sigma_true.2,", Average:",exp(avg_est.2),"\n")
  total_wgt_actual <- sigma_true.1 + sigma_true.2
  total_wgt <- exp(avg_est.1) + exp(avg_est.2)
  cat("Actual ratio (1):",(sigma_true.1/total_wgt_actual),", Average:",(exp(avg_est.1))/total_wgt,"\n")
  cat("Actual ratio (2):",(sigma_true.2/total_wgt_actual),", Average:",(exp(avg_est.2))/total_wgt,"\n")
}

#cat("Matrix normal:\n")
#simulate_matrix_normal(100)

#cat("Matrix-T (1):\n")
#simulate_matT_one(100)

cat("Matrix-T (2):\n")
simulate_matT_two(20)
