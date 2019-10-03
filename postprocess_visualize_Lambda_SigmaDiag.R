source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

#plot_extreme_Lambda(coordinate="mean_x1", level=level, no_indiv=10)
plot_extreme_Lambda(coordinate="mean_x2", level=level, no_indiv=10)

#plot_diag_Sigma(level=level, lrtransform="clr")
