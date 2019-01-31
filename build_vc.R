source("include.R")
suppressMessages(library(driver))
suppressMessages(library(psych))
suppressMessages(library(LaplacesDemon))
suppressMessages(library(mvtnorm))
suppressMessages(library(MCMCpack))
suppressMessages(library(Rcpp))
suppressMessages(library(RcppDE))
suppressMessages(library(RcppEigen))

sourceCpp("dens_optim.cpp")

plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  ggsave(paste(filename,".pdf",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# ============================================================================================
# BUILD VARIANCE COMPONENTS
# ============================================================================================

if(!exists("ilr_data")) {
  load("glom_data_genus.RData")
  f2 <- filter_counts(glom_data, 3, 0.9)

  ilr_data <- NULL
  week_kernel <- NULL
  season_vector <- NULL
  age_vector <- NULL
  group_vector <- NULL
  indiv_vector <- NULL

  #baboons <- c("AMA", "AMO", "ORI")
  #baboons <- c("AMA", "AMO", "BUC", "CHE", "DAG", "EGO", "HOK", "KIW", "NOO", "ORI")
  snames <- unique(read_metadata(f2)$sname)
  #baboons <- snames
  baboons <- snames[sample(length(snames))[1:100]]
  for(b in 1:length(baboons)) {
    cat("Building kernels for baboon",baboons[b],"\n")
    baboon_counts <- subset_samples(f2, sname==baboons[b])
    baboon_metadata <- read_metadata(baboon_counts)
    baboon_log_ratios <- apply(otu_table(baboon_counts)@.Data+0.65, 1, ilr)
    if(is.null(ilr_data)) {
      ilr_data <- baboon_log_ratios
    } else {
      ilr_data <- cbind(ilr_data, baboon_log_ratios)
    }

    nt <- ntaxa(baboon_counts)
    ns <- nsamples(baboon_counts)

    rho <- 1/2
    subset_week_kernel <- matrix(0, nrow=ns, ncol=ns)
    for(i in 1:ns) {
      for(j in 1:ns) {
        if(i == j) {
          subset_week_kernel[i,j] <- 1
        } else {
          d1 <- as.Date(baboon_metadata$collection_date[i])
          d2 <- as.Date(baboon_metadata$collection_date[j])
          # get number of 'weeks difference'
          week_d <- abs(round((d2-d1)/7))[[1]] + 1
          # distance of 1 is < 7 days
          # distance of 2 is < 14 days, etc.
          subset_week_kernel[i,j] <- rho^week_d
        }
      }
    }

    if(is.null(week_kernel)) {
      week_kernel <- subset_week_kernel
    } else {
      prev_r <- dim(week_kernel)[1]
      prev_c <- dim(week_kernel)[2]
      week_kernel <- cbind(week_kernel, matrix(0, nrow=dim(week_kernel)[1], ncol=ns))
      week_kernel <- rbind(week_kernel, matrix(0, nrow=ns, ncol=dim(week_kernel)[2]))
      new_r <- dim(week_kernel)[1]
      new_c <- dim(week_kernel)[2]
      week_kernel[(prev_r+1):new_r,(prev_c+1):new_c] <- subset_week_kernel
    }

    subset_season_vector <- numeric(ns)
    subset_age_vector <- numeric(ns)
    subset_group_vector <- numeric(ns)
    subset_indiv_vector <- numeric(ns)
    for(i in 1:ns) {
      if(baboon_metadata$season[i] == "Wet") {
        subset_season_vector[i] <- "W"
      } else {
        subset_season_vector[i] <- "D"
      }
      if(baboon_metadata$age[i] < 4.5) {
        subset_age_vector[i] <- "J"
      } else if(baboon_metadata$age[i] > 19) {
        subset_age_vector[i] <- "G"
      } else {
        subset_age_vector[i] <- "A"
      }
      subset_group_vector[i] <- baboon_metadata$grp[i]
      subset_indiv_vector[i] <- baboon_metadata$sname[i]
    }
    if(is.null(season_vector)) {
      season_vector <- subset_season_vector
    } else {
      season_vector <- c(season_vector, subset_season_vector)
    }
    if(is.null(age_vector)) {
      age_vector <- subset_age_vector
    } else {
      age_vector <- c(age_vector, subset_age_vector)
    }
    if(is.null(group_vector)) {
      group_vector <- subset_group_vector
    } else {
      group_vector <- c(group_vector, subset_group_vector)
    }
    if(is.null(indiv_vector)) {
      indiv_vector <- subset_indiv_vector
    } else {
      indiv_vector <- c(indiv_vector, subset_indiv_vector)
    }
  }

  #plot_cov(week_kernel, "date_correlation")

  ns_all <- length(season_vector)

  season_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        season_kernel[i,j] <- 1
      } else if(season_vector[i] == "W" && season_vector[j] == "W") {
        season_kernel[i,j] <- 0.1
      } else if(season_vector[i] == "D" && season_vector[j] == "D") {
        season_kernel[i,j] <- 0.2
      } else {
        season_kernel[i,j] <- 0.0
      }
    }
  }
  #plot_cov(season_kernel, "season_correlation")

  age_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        age_kernel[i,j] <- 1
      } else if(age_vector[i] == "J" || age_vector[j] == "J") {
        age_kernel[i,j] <- 0
      } else if(age_vector[i] == "G" && age_vector[j] == "G") {
        age_kernel[i,j] <- 0.2
      } else {
        age_kernel[i,j] <- 0.1
      }
    }
  }
  #plot_cov(age_kernel, "age_correlation")

  group_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        group_kernel[i,j] <- 1
      } else if(group_vector[i] == group_vector[j]) {
        group_kernel[i,j] <- 0.2
      } else {
        group_kernel[i,j] <- 0
      }
    }
  }
  #plot_cov(group_kernel, "group_correlation")

  indiv_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        indiv_kernel[i,j] <- 1
      } else if(indiv_vector[i] == indiv_vector[j]) {
        indiv_kernel[i,j] <- 0.2
      } else {
        indiv_kernel[i,j] <- 0
      }
    }
  }
  #plot_cov(indiv_kernel, "indiv_correlation")
}

# ============================================================================================
# OPTIMIZE VARIANCE COMPONENTS TO TEST DATA
# ============================================================================================

# marginal matrix-T (two component)
logd_mat_t <- function(s) {
  N <- dim(ilr_data)[2]
  P <- dim(ilr_data)[1]
  upsilon <- P + 10
  K <- diag(P)
  A <- exp(s[1])*week_kernel + exp(s[2])*season_kernel + exp(s[3])*group_kernel +
       exp(s[4])*age_kernel + exp(s[5])*indiv_kernel
  logdetA <- 2*sum(log(diag(chol(A))))
  S <- diag(P) + (1/upsilon)*solve(K)%*%(ilr_data)%*%solve(A)%*%t(ilr_data)
  logdetS <- 2*sum(log(diag(chol(S))))
  d <- (P/2)*logdetA + ((upsilon+N+P-1)/2)*logdetS
  return(d)
}

ALT_logd_mat_t <- function(s) {
  return(logd_matrixt(s[1], s[2], s[3], s[4], s[5], ilr_data, week_kernel, season_kernel, group_kernel, age_kernel, indiv_kernel))
}

it <- 1
estimates <- matrix(0, 5, it)
for(i in 1:it) {
  cat("Optimization iteration",i,"\n")
  # out 3x faster if we have R optim call a compiled C++ function
  # could probably speed up further using DEoptim instead of optim here
  res <- optim(par=runif(5), fn=ALT_logd_mat_t, method="L-BFGS-B", lower=rep(0.00001,5), upper=rep(10,5))
  estimates[,i] <- res$par
}
save(estimates, file="Roptim_estimates.RData")

kernels <- c("weekly","seasonal","group","age","individual")
for(i in 1:it) {
  rank <- order(estimates[,i], decreasing=TRUE)
  cat("Iteration 1:")
  for(j in 1:5) {
    cat(" ",kernels[rank[j]]," (",round(estimates[rank[j],i],digits=3),")",sep="")
  }
  cat("\n")
}

# INTERESTING QUESTIONS
# (1) how well is the ordering preserved?
# (2) what proportion of total variance does this model capture? what does the residual variance look like?
# (3) how sensitive is this analysis to choice of VC parameters?
