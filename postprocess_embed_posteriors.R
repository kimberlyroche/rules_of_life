source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]

embed_posteriors(level, which_measure)
