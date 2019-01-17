source("include.R")

filtered <- filter_data()

visualize_groupwise_covariance(filtered, "season", sample=100)
visualize_groupwise_covariance(filtered, "sex", sample=100)
visualize_groupwise_covariance(filtered, "grp", sample=100)
visualize_groupwise_covariance(filtered, "matgrp", sample=100)
visualize_groupwise_covariance(filtered, "sname", sample=100)
visualize_groupwise_covariance(filtered, "age", sample=100)

#perform_mult_timecourse(filtered, c("ACA", "DUI", "CAI", "COB", "DAS"))

#plot_autocorrelation(filtered)
