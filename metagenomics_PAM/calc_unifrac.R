# this utility file can be run as a standalone to calculate (time-consuming) weighted UniFrac distance
# over a set of samples

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

data <- load_and_filter(level=level)
cat("Calculating weighted UniFrac distance at level:",level,"\n")
unifrac_dist <- distance(data, method="unifrac", type="samples", parallel=TRUE)
saveRDS(file.path(relative_path,output_dir,paste0("unifrac_",level,".rds")))
