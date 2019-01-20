source("include.R")

filtered <- filter_data()

visualize_groupwise_covariance(filtered, "sex", sample=1000)
visualize_groupwise_covariance(filtered, "season", sample=1000)
visualize_groupwise_covariance(filtered, "grp", sample=100)
visualize_groupwise_covariance(filtered, "matgrp", sample=100)
visualize_groupwise_covariance(filtered, "sname", sample=100)
visualize_groupwise_covariance(filtered, "age", sample=500)

#perform_mult_timecourse(filtered, c("ACA", "DUI", "CAI", "COB", "DAS"))

#plot_autocorrelation(filtered)
