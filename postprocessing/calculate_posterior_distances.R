# this utility file can be used to (re)calculate posterior distances over either models
# fitted with full posteriors or models fitted to give MAP estimates

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda) (opt: MAP=T/F)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
MAP <- FALSE
if(length(args) == 3) {
  MAP <- as.logical(args[3])
}

# these will be written out to a file, so the assignment is just to prevent output
temp <- calc_posterior_distances(level, which_measure=which_measure, MAP=MAP)
