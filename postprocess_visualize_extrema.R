source("include/R/GP.R")

args <- c("family")
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

plot_extreme_Lambda(coordinate="mean_x1", level=level, no_indiv=1)

plot_extreme_Sigma(coordinate="mean_x1", level=level, no_indiv=1)

plot_diag_Sigma(level=level, lrtransform="clr")

