# this file uses the already-calculated embedding to identify extreme individual with respect
# to the largest principle coordinates of the embedding and to visualize (as heatmaps) either
# the approximate baseline or the approximate average dynamics for these extreme hosts
#
# note: this is probably broken because I've changed the labeling of the embedding coordinates
#       from "mean_x1"; update this

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

# no_indiv=K: the number of extreme individuals to visualize (min=1)
plot_extreme_Lambda(coordinate="mean_x1", level=level, no_indiv=10)

# lrtransform: representation (clr/alr/ilr) to use for the covariance matrix
plot_diag_Sigma(level=level, lrtransform="clr")
