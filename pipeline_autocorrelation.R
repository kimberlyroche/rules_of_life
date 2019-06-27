# plot 16S autocorrelation
# usage: pipeline_autocorrelation.R 16S family 156 weeks FALSE

source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 5) {
  stop("Argument for data type (16S | metagenomics) (level) (lag) (lag units) (resample) missing!\n")
}

data_type <- args[1] # 16S, metagenomics
level <- args[2] # species, genus, family
lag.max <- as.numeric(args[3]) # e.g. 26
lag.units <- args[4] # weeks, months, seasons
resample <- as.logical(args[5])

# read in 16S
glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
# data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE) # 9 is low-count cohort
alr_ref <- ntaxa(data)
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
                             use_alr=TRUE,
                             alr_ref=alr_ref)
if(resample) {
  plot_bounded_autocorrelation(lags,
                               filename=paste("plots/autocorrelation_",lag.max,lag.units,"_",data_type,"_bounded",sep=""),
                               width=10,
                               height=4)
} else {
  plot_mean_autocorrelation(lags,
                            filename=paste("plots/autocorrelation_",lag.max,lag.units,"_",data_type,sep=""),
                            width=10,
                            height=4)
}
