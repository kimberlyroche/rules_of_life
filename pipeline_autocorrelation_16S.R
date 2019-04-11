# plot 16S autocorrelation
# usage: pipeline_autocorrelation_16S.R 1

source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for task missing!\n")
}

task <- as.numeric(args[1])
# 1 : weekly for 26 weeks
# 2 : weekly for 52 weeks
# 3 : monthly for 12 months
# 4 : monthly for 24 months
# 5 : monthly for 36 months
# 6 : seasonal for 5 years

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)

if(task == 1) {
  lags <- calc_autocorrelation(filtered, lag.max=26, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk_",round(100*t),sep=""))
}

if(task == 2) {
  lags <- calc_autocorrelation(filtered, lag.max=52, resample=TRUE, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_52wk")
}

if(task == 3) {
  lags <- calc_autocorrelation(filtered, lag.max=12, resample=TRUE, resample_rate=0.2, date_diff_units="months")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_12mo")
}

if(task == 4) {
  lags <- calc_autocorrelation(filtered, lag.max=24, resample=TRUE, resample_rate=0.2, date_diff_units="months")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_24mo")
}

if(task == 5) {
  lags <- calc_autocorrelation(filtered, lag.max=36, resample=TRUE, resample_rate=0.2, date_diff_units="months")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_36mo")
}

if(task == 6) {
  lags <- calc_autocorrelation(filtered, lag.max=11, resample=TRUE, resample_rate=0.2, date_diff_units="seasons")
  plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_11season")
}
