source("include.R")
library(driver)
library(psych)
library(LaplacesDemon)

plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  ggsave(paste("plots/",filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# pull sample data for two individuals (DUI, ACA)

load("sample_indiv_counts.RData")
md_DUI <- read_metadata(DUI_counts)
md_ACA <- read_metadata(ACA_counts)
ilr_DUI <- apply(otu_table(DUI_counts)@.Data+0.65, 1, ilr)
ilr_ACA <- apply(otu_table(ACA_counts)@.Data+0.65, 1, ilr)
ilr_both <- cbind(ilr_DUI, ilr_ACA)

nt <- dim(ilr_both)[1]
ns_DUI <- dim(ilr_DUI)[2]
ns_ACA <- dim(ilr_ACA)[2]
ns <- ns_DUI + ns_ACA
rm(ilr_DUI)
rm(DUI_counts)
rm(ilr_ACA)
rm(ACA_counts)

# build time kernel (distance per day)

rho <- 1/3
# need "distance" between samples
date_distance <- matrix(0, nrow=ns_DUI+ns_ACA, ncol=ns_DUI+ns_ACA)
for(i in 1:ns_DUI) {
  for(j in 1:ns_DUI) {
    if(i == j) {
      date_distance[i,j] <- 1
    } else {
      dd <- abs(as.Date(md_DUI$collection_date[i]) - as.Date(md_DUI$collection_date[j]))
      if(dd == 0) {
        dd <- 1 # better fix for same-day sampling needed
      }
      dd <- dd[[1]]/50
      date_distance[i,j] <- rho^dd
    }
  }
}
for(i in 1:ns_ACA) {
  for(j in 1:ns_ACA) {
    if(i == j) {
      date_distance[i+ns_DUI,j+ns_DUI] <- 1
    } else {
      dd <- abs(as.Date(md_ACA$collection_date[i]) - as.Date(md_ACA$collection_date[j]))
      if(dd == 0) {
        dd <- 1
      }
      dd <- dd[[1]]/50
      date_distance[i+ns_DUI,j+ns_DUI] <- rho^dd
    }
  }
}
plot_cov(date_distance, "date_correlation_DUI-ACA")

# build individual kernel

indiv_cov <- 0.1
indiv_distance <- matrix(0, nrow=ns_DUI+ns_ACA, ncol=ns_DUI+ns_ACA)
indiv_distance[1:ns_DUI,1:ns_DUI] <- diag(ns_DUI)*(1-indiv_cov) + indiv_cov
indiv_distance[(ns_DUI+1):ns,(ns_DUI+1):ns] <- diag(ns_ACA)*(1-indiv_cov) + indiv_cov
plot_cov(indiv_distance, "indiv_correlation_DUI-ACA")

# build seasonal kernel

season_cov <- 0.1
season_distance <- matrix(0, nrow=ns_DUI+ns_ACA, ncol=ns_DUI+ns_ACA)
for(i in 1:ns_DUI) {
  for(j in 1:ns_DUI) {
    if(i == j) {
      season_distance[i,j] <- 1
    } else {
      if(md_DUI$season[i] == md_DUI$season[j]) {
        season_distance[i,j] <- season_cov
      }
    }
  }
}
for(i in 1:ns_ACA) {
  for(j in 1:ns_ACA) {
    if(i == j) {
      season_distance[i+ns_DUI,j+ns_DUI] <- 1
    } else {
      if(md_ACA$season[i] == md_ACA$season[j]) {
        season_distance[i+ns_DUI,j+ns_DUI] <- season_cov
      }
    }
  }
}
plot_cov(season_distance, "season_correlation_DUI-ACA")

N <- dim(ilr_both)[2]
P <- dim(ilr_both)[1]
mu <- apply(ilr_both, 1, mean)
mu_mat <- matrix(0, nrow=P, ncol=N)
for(i in 1:N) {
  mu_mat[,i] <- mu
}
eta <- ilr_both - mu_mat
image(eta)
K <- cov(t(eta))
image(K)
image(cov(eta))

matT_log_density <- function(s, vc, data, N, P, K) {
  upsilon <- P + 2
  A <- round((exp(s[1])*vc[[1]] + exp(s[2])*vc[[2]] + exp(s[3])*vc[[3]]), digits=10)
#  cat("s1, s2, s3:",s[1],s[2],s[3],"\n")
#  cat("\tlog(det(A)):",log(det(A)),"\n")
#  cat("\tlog(det(solve(A))):",log(det(solve(A))),"\n")
  cat("Tr(A):",tr(A),"\n")
  d <- -(P/2)*log(det(A))-((upsilon+N+P-1)/2)*log(det(diag(P) + solve(K)%*%(data)%*%solve(A)%*%t(data)))
  return(-d)
}

opt_it <- 1
averages <- matrix(0, nrow=3, ncol=opt_it)
for(s in 1:opt_it) {
  cat("Iteration",s,"...\n")
  res <- optim(par=c(runif(1), runif(1), runif(1)),
               vc=list(0.5*date_distance, 0.5*indiv_distance, 0.5*season_distance),
               data=eta, N=N, P=P, K=K, fn=matT_log_density, method="L-BFGS-B")
  total_wgt <- exp(res$par[1]) + exp(res$par[2]) + exp(res$par[3])
  averages[,s] <- c((exp(res$par[1])/total_wgt), (exp(res$par[2])/total_wgt), (exp(res$par[3])/total_wgt))
}
cat("Average percent contribution:",toString(round(apply(averages, 1, mean)*100, 1)),"\n")







