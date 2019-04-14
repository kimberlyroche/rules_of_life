# plot 16S autocorrelation
# usage: pipeline_autocorrelation_16S.R 1

source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 4) {
  stop("Argument for data type (16S | metagenomics) (lag) (lag units) (resample) missing!\n")
}

data_type <- args[1] # 16S, metagenomics
lag.max <- as.numeric(args[2]) # e.g. 26
lag.units <- args[3] # weeks, months, seasons
resample <- as.logical(args[4])

# read in 16S
glom_data <- load_glommed_data(level="species", replicates=TRUE)
data <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(data)

resample_rate = 0.2

if(data_type == "metagenomics") {
  # read in metagenomics/PiPhillin
  data.piphillin <- read_metagenomics(metadata)
  metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)
  data <- data.piphillin
  metadata <- metadata.metagenomics
}

if(resample) {
  lags <- calc_autocorrelation(data, metadata, lag.max=lag.max, date_diff_units=lag.units, resample=resample, resample_rate=0.2)
  plot_bounded_autocorrelation(lags, filename=paste("plots/autocorrelation_",lag.max,lag.units,"_",data_type,"_bounded",sep=""))
} else {
  lags <- calc_autocorrelation(data, metadata, lag.max=lag.max, date_diff_units=lag.units, resample=resample, resample_rate=0.2)
  plot_mean_autocorrelation(lags, filename=paste("plots/autocorrelation_",lag.max,lag.units,"_",data_type,sep=""))
}
