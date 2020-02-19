# this file plots visualizations for a given host in the form of bar/composition plots, fits for
# specified log relative abundances of microbes, and heatmaps of approximate average correlation between
# log relative abundances of microbes

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))
source(file.path(relative_path,"include/R/visualization.R"))
source(file.path(relative_path,"include/R/data_transform.R"))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (coordinate) (sname1) (sname2) ...", call.=FALSE)
}
level <- args[1]
predict_coords <- c(as.numeric(args[2]))
individuals <- c()
if(length(args) > 2) {
  for(i in 3:length(args)) {
    individuals <- c(individuals, args[i])
  }
}
if(length(individuals) == 0) {
  individuals <- over_40[sample(1:length(over_40))[1:10]]
}

# timecourse=T/F: plot bar/composition plots
# covcor=T/F: plot a heatmap of this individuals' mean posterior covariance and correlation over taxa
# predict_coords: we can pass (predict_coords=predict_coords) or not pass (predict_coords=NULL) log ratio
#                 coordinates to indicate to the function we want to plot Gaussian process predictions 
#                 for these log relative abundances

plot_ribbons_individuals(individuals, level, timecourse=TRUE, covcor=FALSE, predict_coords=NULL)
