level <- "family"

check_indivs <- c("ECH", "LAV", "LOB", "LYE", "MON", "NOB", "ONY", "OPH", "QUA", "TAL", "VET", "YAI")
#check_indivs <- NULL

if(is.null(check_indivs)) {
  pattern_str <- "*_bassetfit.rds"
  regexpr_str <- "_bassetfit.rds"
  fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern=pattern_str, full.names=TRUE, recursive=FALSE)
  for(fm in fitted_models) {
    fit <- readRDS(fm)$fit
    if(is.nan(fit$logMarginalLikelihood)) {
      cat(fm,"\n")
      cat("\tLL:",fit$logMarginalLikelihood,"\n\n")
    }
  }
} else {
  for(indiv in check_indivs) {
    fm <- paste0("subsetted_indiv_data/",level,"/",indiv,"_bassetfit.rds")
    fit <- readRDS(fm)$fit
    cat(fm,"\n")
    cat("\tLL:",fit$logMarginalLikelihood,"\n\n")
  }
}
