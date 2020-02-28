# the file independently fits the model to the individuals named in the sname_list variable ultimately included
# from `include/R/general.R` (through `GP.R`)
#
# output is a stray::pibble_fit object with posterior samples Eta, Lambda, and Sigma (covariance over taxa);
# these are saved in the location specified by the model_dir variable again in `include/R/general.R`

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

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
  se_sigma <- as.numeric(args[4])
  per_sigma <- as.numeric(args[5])
  wn_sigma <- as.numeric(args[6])
} else {
  # defaults; these were determined by using `formalize_parameters.R` in the preprocessing folder
  # these are heuristic values based on filtering 5-counts, at 20% min sample representation
  #
  # the ALR reference is taken to be the taxon with median log abundance, just to prevent it being something very rare
  #
  # wn_sigma species the scale of the (unused) white noie kernel
  if(level == "phylum") {
    se_sigma <- sqrt(2.052)
    per_sigma <- sqrt(0.228)
    wn_sigma <- sqrt(0)
    alr_ref <- 9
  }
  if(level == "class") {
    se_sigma <- sqrt(2.177)
    per_sigma <- sqrt(0.242)
    wn_sigma <- sqrt(0)
    alr_ref <- 15
  }
  if(level == "order") {
    se_sigma <- sqrt(2.52)
    per_sigma <- sqrt(0.28)
    wn_sigma <- sqrt(0)
    alr_ref <- 17
  }
  if(level == "family") {
    se_sigma <- sqrt(2.586)
    per_sigma <- sqrt(0.287)
    wn_sigma <- sqrt(0)
    alr_ref <- 22
  }
  if(level == "genus") {
    wn_sigma <- 0
    se_sigma <- sqrt(2.632)
    per_sigma <- sqrt(0.292)
    alr_ref <- 59
    # alt filter: 5-count x 50%
    # se_sigma <- sqrt(2.601)
    # per_sigma <- sqrt(0.289)
    # alr_ref <- 28
    # alt filter: 20-count x 20%
    # se_sigma <- sqrt(2.677)
    # per_sigma <- sqrt(0.297)
    # alr_ref <- 47
    # alt filter: 20-count x 50%
    # se_sigma <- sqrt(2.273)
    # per_sigma <- sqrt(0.253)
    # alr_ref <- 24
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

fit_GP(baboon, level, se_sigma=se_sigma, per_sigma=per_sigma, wn_sigma=wn_sigma,
       dd_se=dd_se, save_append=save_append, date_lower_limit=NULL, date_upper_limit=NULL,
       alr_ref=alr_ref, mean_only=mean_only,
       max_iter=20000, eps_f=1e-11, eps_g=1e-5)
