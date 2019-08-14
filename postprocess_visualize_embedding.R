source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda) (label type)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
label_type <- args[3]

plot_ordination(level, which_measure, label_type, legend=TRUE)
