source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop(paste("Arguments: (sname) (level) (OPT: SE scale) (OPT: PER scale) ",
            "(OPT: days to 0.1 correlation decay) (OPT: save filename slug)"), call.=FALSE)
}
baboon <- args[1]
level <- args[2]
if(length(args) >= 3) {
  se_weight <- as.numeric(args[3])
  per_weight <- se_weight*as.numeric(args[4])
} else {
  # defaults
  se_weight <- 2
  per_weight <- 0.25
}
if(length(args) >= 5) {
  dd_se <- as.numeric(args[5])
} else {
  dd_se <- 90
}
if(length(args) >= 6) {
  save_append <- paste0("_",args[6])
} else {
  save_append <- ""
}

fit_GP(baboon, level, se_weight=se_weight, per_weight=per_weight,
       dd_se=dd_se, save_append=save_append, date_lower_limit=NULL, date_upper_limit=NULL)
