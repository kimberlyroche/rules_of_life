library(mvtnorm)
library(driver)
library(ggplot2)
library(stringr)

source("include/R/visualization.R")

baboons <- c("VET", "GAN")
for(b in baboons) {
  cat("Simulating individual",b,"...\n")
  # get a baseline
  #temp <- readRDS(paste0("output/model_fits/genus/LEI_bassetfit.rds"))
  temp <- readRDS(paste0("output/model_fits/genus/",b,"_bassetfit.rds"))
  baseline <- t(temp$alr_ys)
  baseline <- apply(baseline, 1, mean)
  # get "dynamics"
  #temp <- readRDS(paste0("output/model_fits/genus/",b,"_bassetfit.rds"))
  Sigma <- cov2cor(apply(temp$fit$Sigma, c(1,2), mean))
  p <- nrow(Sigma)
  n <- 100
  simulation <- matrix(0, p, n)
  simulation[,1] <- baseline
  for(i in 2:n) {
    simulation[,i] <- rmvnorm(1, mean=simulation[,(i-1)], sigma=Sigma)
  }
  prop <- t(driver::alrInv(t(simulation)))

  df <- as.data.frame(prop)
  colnames(df) <- paste("Sample", str_pad(1:n, 3, pad="0"))
  rownames(df) <- paste("OTU", str_pad(1:(p+1), 3, pad="0"))
  df2 <- df %>% gather("sample", "value")
  df2$OTU <- rep(seq(1:(p+1)), n)
  df2$OTU <- as.factor(df2$OTU)
  categories <- unique(df2$OTU)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(categories))
  pl <- ggplot(df2, aes(x=sample, y=value, fill=OTU)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8)) +
    theme(legend.position="none")
  ggsave(paste0(b,"_simulation.png"), plot=pl, scale=1.5, dpi=100, width=10, height=6, units="in")
}
