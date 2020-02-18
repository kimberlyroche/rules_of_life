library(phyloseq)

source("include/R/general.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

data <- load_and_filter(level=level)
cat("Calculating weighted UniFrac distance at level:",level,"\n")
unifrac_dist <- distance(data, method="unifrac", type="samples", parallel=TRUE)
saveRDS(paste0(output_dir,"unifrac_",level,".rds"))
