library(Rcpp)
library(ggplot2)

source("include.R")
sourceCpp("mat_dist.cpp")

filenames <- c("test_ALR", "test_ILR", "test_CLR")
for(fn in filenames) {
  Sigma <- readRDS(paste0(fn,".rds"))$fit$Sigma
  if(fn == "test_CLR") {
    truncSigma <- array(NA, dim=c(dim(fit_obj$fit$Sigma)[1]-1,dim(fit_obj$fit$Sigma)[2]-1, dim(fit_obj$fit$Sigma)[3]))
    for(i in 1:dim(fit_obj$fit$Sigma)[3]) {
      truncSigma[,,i] <- fit_obj$fit$Sigma[1:(dim(fit_obj$fit$Sigma)[1]-1),1:(dim(fit_obj$fit$Sigma)[2]-1),i]
    }
    Sigma <- truncSigma
  }
  dim(Sigma) <- c(dim(Sigma)[1], dim(Sigma)[2]*dim(Sigma)[3])
  cat("Calculating distances...\n")
  d <- mat_dist(Sigma, 1, 100)
  cat("Embedding...\n")
  fit <- cmdscale(d, eig=TRUE, k=2)
  df <- data.frame(x=fit$points[,1], y=fit$points[,2], labels=as.factor(1:100))
  if(fn == "test_CLR") {
    df$x <- df$x
    df$y <- -df$y
  }
  p <- ggplot() +
         geom_text(data=df, aes(x=x, y=y, label=labels), color="black", fontface="bold")
  ggsave(paste0(fn,".png"), plot=p, scale=1.5, width=6, height=6, units="in", dpi=72)
}
