source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda) (opt: MAP=T/F) (opt: sname)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
indiv <- NULL
MAP <- FALSE
if(length(args) == 3) {
  MAP <- as.logical(args[3])
}
if(length(args) == 4) {
  indiv <- args[4]
}

embed_posteriors(level, which_measure, MAP=MAP, indiv=indiv)

# Euclidean distance in ILR
# embed_posteriors_alt(level, indiv=indiv)
