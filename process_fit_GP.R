source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop(paste("Arguments: (sname) (level) (OPT: SE scale) (OPT: PER scale) (OPT: noise scale) ",
            "(OPT: days to 0.1 correlation decay) (OPT: save filename slug)"), call.=FALSE)
}
baboon <- args[1]
level <- args[2]
if(length(args) >= 3) {
  se_weight <- as.numeric(args[3])
  per_weight <- as.numeric(args[4])
  wn_weight <- as.numeric(args[5])
} else {
  # defaults
  if(level == "phylum") {
    se_weight <- sqrt(1.956)
    per_weight <- sqrt(0.217)
    wn_weight <- sqrt(0.233)
  }
  if(level == "family") {
    se_weight <- sqrt(2.090)
    per_weight <- sqrt(0.232)
    wn_weight <- sqrt(0.381)
  }
  if(level == "genus") {
    wn_weight <- sqrt(0.354)
    se_weight <- sqrt(2.181)
    per_weight <- sqrt(0.242)
  }
}
if(length(args) >= 6) {
  dd_se <- as.numeric(args[6])
} else {
  dd_se <- 90
}
if(length(args) >= 7) {
  save_append <- paste0("_",args[7])
} else {
  save_append <- ""
}

fit_GP(baboon, level, se_weight=se_weight, per_weight=per_weight, wn_weight=wn_weight,
       dd_se=dd_se, save_append=save_append, date_lower_limit=NULL, date_upper_limit=NULL,
       max_iter=20000, eps_f=1e-11, eps_g=1e-5)
