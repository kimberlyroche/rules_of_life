source("include/R/GP.R")
source("include/R/visualization.R")
source("include/R/data_transform.R")

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

# update plot_ribbons_individual to facet_wrap multiple individuals!

plot_ribbons_individuals(individuals, level, timecourse=FALSE, covcor=FALSE, predict_coords=predict_coords)

#plot_ribbons_individuals(individuals, level, timecourse=FALSE, covcor=TRUE, predict_coords=NULL)

#plot_ribbons_individuals(individuals, level, timecourse=TRUE, covcor=TRUE, predict_coords=NULL)
