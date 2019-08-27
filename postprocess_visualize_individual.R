source("include/R/GP.R")
source("include/R/visualization.R")
source("include/R/data_transform.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (sname1) (sname2) ...", call.=FALSE)
}
level <- args[1]
individuals <- c()
for(i in 2:length(args)) {
  individuals <- c(individuals, args[i])
}

plot_ribbons_individuals(individuals, level, timecourse=TRUE, covcor=FALSE, predict_coords=c(1,10))
