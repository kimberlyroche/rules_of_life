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

# read in 16S
glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(filtered)

# read in metagenomics/PiPhillin
data.piphillin <- read_metagenomics(metadata)

if(task == 1) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=26, resample=F, resample_rate=0.2)
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=26, resample=T, resample_rate=0.2)
  #plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk",sep=""))
  
  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=26, resample=F, resample_rate=0.2)
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_26wk_metagenomics",sep=""))
}

if(task == 2) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=52, resample=F, resample_rate=0.2)
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_52wk",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=52, resample=T, resample_rate=0.2)
  #plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_52wk")
  
  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=52, resample=F, resample_rate=0.2)
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_52wk_metagenomics",sep=""))
}

if(task == 3) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=12, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_12mo",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=12, resample=T, resample_rate=0.2, date_diff_units="months")
  #plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_12mo")

  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=12, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_12mo_metagenomics",sep=""))
}

if(task == 4) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=24, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_24mo",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=24, resample=T, resample_rate=0.2, date_diff_units="months")
  #plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_24mo")
  
  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=24, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_24mo_metagenomics",sep=""))
}

if(task == 5) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=36, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_36mo",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=36, resample=T, resample_rate=0.2, date_diff_units="months")
  #plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_36mo")
  
  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=36, resample=F, resample_rate=0.2, date_diff_units="months")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_36mo_metagenomics",sep=""))
}

if(task == 6) {
  lags <- calc_autocorrelation(filtered, metadata, lag.max=11, resample=F, resample_rate=0.2, date_diff_units="seasons")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_11season",sep=""))
  #lags <- calc_autocorrelation(filtered, metadata, lag.max=11, resample=T, resample_rate=0.2, date_diff_units="seasons")
  #plot_bounded_autocorrelation(lags, filename="plots/autocorrelation_11season")

  lags <- calc_autocorrelation(data.piphillin, metadata, lag.max=11, resample=F, resample_rate=0.2, date_diff_units="seasons")
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_11season_metagenomics",sep=""))
}








