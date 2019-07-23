# plot time course and mean covariance (ILR) for a given individual

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript 02_mean_Sigma.R ACA family", call.=FALSE)
}
baboon <- args[1]
level <- args[2]

cat("Baboon:",baboon,", level:",level,"\n")

source("include.R")

glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=10, sample_threshold=0.66, verbose=TRUE)
indiv_data <- subset_samples(data, sname==baboon)
cat("\tPlotting timecourse...\n")
# this function already prepends with 'plot' -- watch out!
plot_timecourse_phyloseq(indiv_data, paste0("basset/",level,"/",baboon,"_timecourse"), gapped=FALSE,
                                   legend=TRUE, legend_level=level)
cat("\tPlotting mean covariance/correlation...\n")
# load Sigma samples
fn <- paste0("subsetted_indiv_data/",level,"/",baboon,"_bassetfit.RData")
load(fn)
cat("\t\tMean trace:",mean(apply(Sigma, 3, function(x) sum(diag(x)))),"\n")
temp <- apply(Sigma, c(1,2), mean)
png(paste0("plots/basset/",level,"/",baboon,"_mean_cov.png"))
image(temp)
dev.off()
png(paste0("plots/basset/",level,"/",baboon,"_mean_corr.png"))
image(cov2cor(temp))
dev.off()
