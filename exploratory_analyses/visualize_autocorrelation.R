# this file produces autocorrelation plots at specified lag, with or without bootstrapping
# of samples to approximate uncertainty

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 5) {
  stop("Arguments: (16S | metagenomics) (level) (lag) (lag units) (resample)\n")
}
data_type <- args[1] # 16S, metagenomics
level <- args[2] # species, genus, family
lag.max <- as.numeric(args[3]) # e.g. 26
lag.units <- args[4] # weeks, months, seasons
resample <- as.logical(args[5])

data <- load_and_filter(level)
metadata <- read_metadata(data)

if(data_type == "metagenomics") {
  # read in metagenomics/PiPhillin
  data.piphillin <- read_metagenomics(metadata)
  metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)
  data <- data.piphillin
  metadata <- metadata.metagenomics
}

lags <- calc_autocorrelation(data,
                             metadata,
                             lag.max=lag.max,
                             date_diff_units=lag.units,
                             resample=resample,
                             use_lr="alr")
if(resample) {
  plot_bounded_autocorrelation(lags,
                               filename=file.path(relative_path,plot_dir,paste0("autocorrelation_",lag.max,lag.units,"_",data_type,"_bounded")),
                               width=10,
                               height=4)
} else {
  plot_mean_autocorrelation(lags,
                            filename=file.path(relative_path,plot_dir,paste0("autocorrelation_",lag.max,lag.units,"_",data_type)),
                            width=10,
                            height=4)
}
