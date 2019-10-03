source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda) (opt: sname)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
indiv <- NULL
if(length(args) == 3) {
  indiv <- args[3]
}

embed_posteriors(level, which_measure, indiv=indiv)
