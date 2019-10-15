library(stray)
library(ggplot2)

source("include/R/general.R")

level <- "genus"

#model_count_limit <- 5

# read all models
#fns <- list.files(path=paste0(model_dir, level))
fns <- list.files(path="output/model_fits/genus_prev")
models <- list()
ref_model <- NULL
for(i in 1:length(fns)) {
#for(i in 1:model_count_limit) {
  cat("Reading model:",fns[i],"\n")
  #temp <- readRDS(paste0(model_dir, level, "/", fns[i]))
  temp <- readRDS(paste0("output/model_fits/genus_prev/", fns[i]))
  if(is.null(ref_model)) {
    ref_model <- temp
  }
  temp.fit <- to_clr(temp$fit)
  temp.Sigma <- temp.fit$Sigma
  for(k in 1:dim(temp.Sigma)[3]) {
    temp.Sigma[,,k] <- cov2cor(temp.Sigma[,,k])
  }
  models[[i]] <- apply(temp.Sigma, c(1,2), mean)
}
cat("Models read.\n")

cd <- data.frame(m1=c(), m2=c(), microbe=c(), cor=c(), avg_log_abundance=c())
model_pairs <- combn(length(models), 2)
for(i in 1:ref_model$fit$D) {
  cat("Computing correlation over dynamics for microbe",i,"\n")
  for(j in 1:length(models)) {
    for(k in 1:length(models)) {
      if(j != k) {
        cd <- rbind(cd,
              data.frame(m1=j,
                         m2=k,
                         microbe=i,
                         cor=cor(models[[j]][i,], models[[k]][i,]),
                         avg_log_abundance=mean(log(ref_model$Y[i,] + 0.5))))
      }
    }
  }
}
saveRDS(cd,"tax_dynamics_corr.rds")

p <- ggplot(cd, aes(cor, color=avg_log_abundance)) +
  geom_density() +
  facet_wrap(~ microbe) +
  scale_color_gradient(low="blue", high="red") +
  labs(fill = "Mean log(abundance)")
ggsave("tax_relationship_corr_distro.png", dpi=100, scale=1.5, units="in", width=15, height=9)

