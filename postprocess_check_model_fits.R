# check stray::basset model fits
# occasionally the model converges to an extremely large scale for Sigma
# need to debug; for now, just identify the bad fits and re-run
# these fits have vanishingly low likelihood (NaN), so search for that explicitly, either
#   just within the first for individuals identified in check_indivs or for all fitted
#   models if check_indivs is NULL

source("include/R/GP.R")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  cat("Arguments: (level)")
  quit()
}

level <- args[1]

#check_indivs <- c("ECH", "LAV", "LOB", "LYE", "MON", "NOB", "ONY", "OPH", "QUA", "TAL", "VET", "YAI")
check_indivs <- NULL

if(is.null(check_indivs)) {
  fitted_models <- get_fitted_modellist(level)$fitted_models
  for(fm in fitted_models) {
    fit <- readRDS(fm)$fit
    if(is.nan(fit$logMarginalLikelihood)) {
      cat(fm,"\n")
      cat("\tLL:",fit$logMarginalLikelihood,"\n\n")
    }
  }
} else {
  for(indiv in check_indivs) {
    fm <- paste0(model_dir,level,"/",indiv,"_bassetfit.rds")
    fit <- readRDS(fm)$fit
    cat(fm,"\n")
    cat("\tLL:",fit$logMarginalLikelihood,"\n\n")
  }
}
