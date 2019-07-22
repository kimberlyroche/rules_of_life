args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript 02_mean_Sigma.R family 90_1", call.=FALSE)
}
level <- args[1]
if(length(args) > 1) {
  tag <- args[2]
} else {
  tag <- NULL
}

pattern_str <- "*_bassetfit.RData"
regexpr_str <- "_bassetfit.RData"
if(!is.null(tag)) {
  pattern_str <- paste0("*_bassetfit_",tag,".RData")
  regexpr_str <- paste0("_bassetfit_",tag,".RData")
}
fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern=pattern_str, full.names=TRUE, recursive=FALSE)
individuals <- sapply(fitted_models, function(x) { idx <- regexpr(regexpr_str, x); return(substr(x, idx-3, idx-1)) } )
names(individuals) <- NULL

for(i in 1:length(individuals)) {
  cat("Processing baboon",individuals[i],"...\n")
  fn <- paste0("subsetted_indiv_data/",level,"/",individuals[i],regexpr_str)
  load(fn)
  cat("Mean trace:",mean(apply(Sigma, 3, function(x) sum(diag(x)))),"\n")
  cat("Sigma has",dim(Sigma)[3],"samples\n")
  temp <- apply(Sigma, c(1,2), mean)
  png(paste0("plots/basset/",level,"/",individuals[i],"_meanpostcorr.png"))
  image(cov2cor(temp))
  dev.off()
}
