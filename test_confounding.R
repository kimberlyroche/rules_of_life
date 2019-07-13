source("include.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  # argument is model tags
  stop("Usage: Rscript test_confounding.R 180_1", call.=FALSE)
}
tag <- args[1]

level <- "family"
glom_data <- load_glommed_data(level=level, replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.66)

fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern=paste0("*_bassetfit_",tag,".RData"), full.names=TRUE, recursive=FALSE)
individuals <- sapply(fitted_models, function(x) { idx <- regexpr(paste0("_bassetfit_",tag,".RData"), x); return(substr(x, idx-3, idx-1)) } )
cat(individuals)
names(individuals) <- NULL

sample_count <- c()
sample_density <- c()
var_scale <- c()
off_diag_net <- c()
off_diag_mass <- c()
max_d <- -Inf
max_d_indiv <- NULL
max_d_days <- 0
min_d <- Inf
min_d_indiv <- NULL
min_d_days <- 0
for(indiv in individuals) {
  cat("Parsing individual",indiv,"...\n")
  indiv_subset <- subset_samples(filtered, sname==indiv)
  this_count <- phyloseq::nsamples(indiv_subset)
  sample_count <- c(sample_count, this_count)
  metadata <- sample_data(indiv_subset)
  dates_col <- sort(unique(metadata$collection_date))
  col_range <- date_to_num(dates_col[length(dates_col)], dates_col[1])
  this_density <- this_count/col_range
  sample_density <- c(sample_density, this_density)
  if(sample_density < min_d) {
    min_d <- sample_density
    min_d_indiv <- indiv
    min_d_days <- col_range
  }
  if(sample_density > max_d) {
    max_x <- sample_density
    max_d_indiv <- indiv
    max_d_days <- col_range
  }
  load(paste0("subsetted_indiv_data/family/",indiv,"_bassetfit_",tag,".RData"))
  var_scale <- c(var_scale, sum(diag(apply(Sigma, c(1,2), mean))))
  off_diag <- apply(Sigma, 3, function(x) { y=cov2cor(x); diag(y)=0; mean(y) } )
  off_diag_net <- c(off_diag_net, mean(off_diag))
  off_diag_mass <- c(off_diag_mass, mean(abs(off_diag)))
}

cat("Individual with max density:",max_d_indiv," (",max_d_days,")\n")
cat("Individual with min density:",min_d_indiv," (",min_d_days,")\n\n")

#cat("Range of sample counts:",min(sample_count),",",max(sample_count),"\n")
#cat("Range of sample densities:",min(sample_density),",",max(sample_density),"\n\n")

cat("Sample count vs. total variance:",cor(sample_count, var_scale),"\n")
cat("Sample density vs. total variance:",cor(sample_density, var_scale),"\n\n")

cat("Sample count vs. signed interations:",cor(sample_count, off_diag_net),"\n")
cat("Sample density vs. signed interactions:",cor(sample_density, off_diag_net),"\n\n")

cat("Sample count vs. unsigned interations:",cor(sample_count, off_diag_mass),"\n")
cat("Sample density vs. unsigned interactions:",cor(sample_density, off_diag_mass),"\n")
