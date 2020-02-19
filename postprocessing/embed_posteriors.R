# this file pulls the fitted models, calculates distances between posterior samples (if they
# don't already exist) and embeds them (using MDS)
#
# an optional host short name can be provided (variable indiv below); not sure why I included this?

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

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
