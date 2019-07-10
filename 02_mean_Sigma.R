args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript 02_mean_Sigma.R family", call.=FALSE)
}
level <- args[1]

fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern="*_bassetfit_.RData", full.names=TRUE, recursive=FALSE)
individuals <- sapply(fitted_models, function(x) { idx <- regexpr("_bassetfit_.RData", x); return(substr(x, idx-3, idx-1)) } )
names(individuals) <- NULL

for(i in 1:length(individuals)) {
  cat("Processing baboon",individuals[i],"...\n")
  load(paste0("subsetted_indiv_data/",level,"/",individuals[i],"_bassetfit_.RData"))
  cat("Mean trace:",mean(apply(Sigma, 3, function(x) sum(diag(x)))),"\n")
  cat("Sigma has",dim(Sigma)[3],"samples\n")
  temp <- apply(Sigma, c(1,2), mean)
  png(paste0("plots/basset/",level,"/Sigma_mean_",individuals[i],".png"))
  image(temp)
  dev.off()
}
