best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")

for(i in 1:length(best_sampled)) {
  cat("Processing baboon",best_sampled[i],"...\n")
  load(paste0("subsetted_indiv_data/",best_sampled[i],"_bassetfit.RData"))
  cat("Mean trace:",mean(apply(fit_obj$fit$Sigma, 3, function(x) sum(diag(x)))),"\n")
  cat("Sigma has",dim(fit_obj$fit$Sigma)[3],"samples\n")
  temp <- apply(fit_obj$fit$Sigma, c(1,2), mean)
  png(paste0("plots/basset/crude_Sigma_mean_",best_sampled[i],".png"))
  image(temp)
  dev.off()
}
