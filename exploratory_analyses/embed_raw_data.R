# TBD: the idea here is to directly embed the data and see if there is any evidence of "modes" at this level

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

sourceCpp(file.path(relative_path,"include/cpp/uncollapse_simple.cpp"))

level <- "family"
data <- load_and_filter(level)
md <- sample_data(data) # these are individuals with >= 50 (40?) samples

# (1) log relative abundance transform each individual's samples

# (2) detrend individually for time using uncollapse_simple(eta, X, Theta, Gamma, Xi, upsilon)

# (3) calculate empirical covariance

# (4) calculate Riemannian distance

# (5) embed using MDS
  
  