source("include.R")

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("Argument for filtering stringency missing!\n")
}

filter_percent <- as.numeric(args[1])

# not currently in use; need to add it to include.R at some point
# plots crude variance components (per-individual) for diagnostic purposes
plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  #ggsave(paste(filename,".pdf",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
  ggsave(paste(filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# ============================================================================================
# BUILD VARIANCE COMPONENTS
# ============================================================================================

load("glom_data_genus.RData")
filtered <- filter_counts(glom_data, 3, filter_percent)
estimate_variance_components(filtered=filtered, optim_it=1)
