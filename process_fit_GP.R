source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop(paste("Arguments: (sname) (level) (OPT: sample #) (OPT: SE scale) (OPT: PER scale) (OPT: noise scale) ",
            "(OPT: days to 0.1 correlation decay) (OPT: save filename slug)"), call.=FALSE)
}
baboon <- args[1]
level <- args[2]
mean_only <- FALSE
alr_ref <- NULL
if(length(args) >= 3) {
  mean_only <- as.logical(args[3])
}
if(length(args) >= 4) {
  se_weight <- as.numeric(args[4])
  per_weight <- as.numeric(args[5])
  wn_weight <- as.numeric(args[6])
} else {
  # defaults
  if(level == "phylum") {
    se_weight <- sqrt(2.052)
    per_weight <- sqrt(0.228)
    wn_weight <- sqrt(0)
    alr_ref <- 9
  }
  if(level == "family") {
    se_weight <- sqrt(2.586)
    per_weight <- sqrt(0.287)
    wn_weight <- sqrt(0)
    alr_ref <- 22
  }
  if(level == "genus") {
    wn_weight <- 0
    # filter 5, 20%
    se_weight <- sqrt(2.632)
    per_weight <- sqrt(0.292)
    alr_ref <- 59
    # filter 5, 50%
    #se_weight <- sqrt(2.601)
    #per_weight <- sqrt(0.289)
    #alr_ref <- 28
    # filter 20, 20%
    #se_weight <- sqrt(2.677)
    #per_weight <- sqrt(0.297)
    #alr_ref <- 47
    # filter 20, 50%
    #se_weight <- sqrt(2.273)
    #per_weight <- sqrt(0.253)
    #alr_ref <- 24
  }
}

if(length(args) >= 7) {
  dd_se <- as.numeric(args[7])
} else {
  dd_se <- 90
}
if(length(args) >= 8) {
  save_append <- paste0("_",args[8])
} else {
  save_append <- ""
}

fit_GP(baboon, level, se_weight=se_weight, per_weight=per_weight, wn_weight=wn_weight,
       dd_se=dd_se, save_append=save_append, date_lower_limit=NULL, date_upper_limit=NULL,
       alr_ref=alr_ref, mean_only=mean_only,
       max_iter=20000, eps_f=1e-11, eps_g=1e-5)
