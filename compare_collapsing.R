source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for agglomeration-level missing!\n")
}
if(length(args)==1) {
  stop("Argument for filter-level missing!\n")
}

l <- args[1] # "na", "species", "genus"
f <- as.numeric(args[2]) # 0.2, 0.5, 0.9


# assume the existence of:
#   glom_data_genus.RData
#   glom_data_species.RData
#   data/emp_baboon_NewFiltr.RDS

glom_data <- readRDS("data/emp_baboon_NewFiltr.RDS")

original_counts <- sum(otu_table(glom_data)@.Data)
original_approx_zero_prop <- 1 - get_tiny_counts(glom_data, 1, use_ns=1000)
cat("Original total counts:",original_counts,"\n")
cat("Original approx. zero proportion:",original_approx_zero_prop,"\n")

if(l == "species") {
  load("glom_data_species.RData")
}
if(l == "genus") {
  load("glom_data_genus.RData")
}

filtered <- filter_data(sample_threshold=f, data=glom_data)

# (1) compare retained counts
retained_counts <- sum(otu_table(filtered)@.Data)
cat("New total counts:",retained_counts,"\n")

# (2) compare reduction in zero counts
approx_zero_prop <- 1 - get_tiny_counts(filtered, 1, use_ns=1000)
cat("Approx. zero proportion:",approx_zero_prop,"\n")

# (3) compare autocorrelation plot at 52-weeks
lags <- calc_autocorrelation(filtered, lag.max=52, resample=TRUE, resample_rate=0.2)
plot_bounded_autocorrelation(lags, filename=paste("comparing_collapse_params/compare_AC_f",round(100*f),"_l-",l,sep=""))

# (4) compare variance component output
estimate_variance_components(filtered=filtered, optim_it=2)
