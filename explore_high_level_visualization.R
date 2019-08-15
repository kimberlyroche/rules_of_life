source("include/R/general.R")

# NOTE: refactored and untested

# ====================================================================================================================
# get noise (1) from replicates and (2) overall
# ====================================================================================================================

level <- "species"
data <- load_and_filter(level)
replicates <- subset_samples(data, sample_status==2)
metadata <- sample_data(replicates)

rep_sid <- unique(metadata$sid)
var_clr <- c()
var_alr <- c()
for(rs in rep_sid) {
	rep_subset <- subset_samples(replicates, sid==rs)
	counts <- otu_table(rep_subset)@.Data
	sample_counts <- apply(counts, 1, sum)
	if(nrow(counts) > 2 & min(sample_counts) < 5000) {
		exclude <- which.min(sample_counts)
		cat("Excluding",exclude,"with counts =",sample_counts[exclude],"!\n")
		counts <- counts[-exclude,]
	}
	alr_counts <- driver::alr(counts + 0.5)
	clr_counts <- driver::clr(counts + 0.5)
	var_alr <- c(var_alr, apply(alr_counts, 2, var))
	var_clr <- c(var_clr, apply(clr_counts, 2, var))
}
cat("CLR replicate variance:",mean(var_clr),"\n")
cat("ALR replicate variance:",mean(var_alr),"\n")

counts <- otu_table(filtered)@.Data
clr_counts <- driver::clr(counts + 0.5)
alr_counts <- driver::alr(counts + 0.5)
var_clr <- apply(clr_counts, 2, var)
var_alr <- apply(alr_counts, 2, var)
cat("CLR variance:",mean(var_clr),"\n")
cat("ALR variance:",mean(var_alr),"\n")

# ====================================================================================================================
# plot sample overview: all samples with time on x-axis x individual on y-axis
# ====================================================================================================================

# ====================================================================================================================
# read in metagenomics
# ====================================================================================================================

# read metagenomics
data.piphillin <- read_metagenomics(metadata)
metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)

if(FALSE) {

color_metagenomics <- TRUE
  
cat("Earliest sample:", min(metadata$collection_date), "\n")
cat("Latest sample:", max(metadata$collection_date), "\n")

min_samples <- 10 # drop individuals with fewer than this many samples
individuals <- unique(metadata$sname)
samples <- list()
for(i in 1:length(individuals)) {
  indiv_samples <- subset_samples(metadata, (sname %in% c(individuals[i])))
  indiv_collection_dates <- indiv_samples$collection_date
  indiv_sample_id <- indiv_samples$sample_id
  indiv_samples <- as.vector(sapply(indiv_collection_dates, date_to_num))
  if(length(indiv_samples) >= min_samples) {
    samples[[individuals[i]]] <- list(indiv_samples, indiv_sample_id)
  }
}

metagenomics.sample_ids <- metadata.metagenomics[,"sample_id"]$sample_id

plot.data <- data.frame(indiv=NULL, sample=NULL, metagenomics=NULL)
for(ind in names(samples)) {
  if(color_metagenomics) {
    new.data <- data.frame(indiv=rep(ind, length(samples[[ind]][[1]])),
                           sample=samples[[ind]][[1]], metagenomics=as.factor(samples[[ind]][[2]] %in% metagenomics.sample_ids))
  } else {
    new.data <- data.frame(indiv=rep(ind, length(samples[[ind]][[1]])),
                           sample=samples[[ind]][[1]])
  }
  plot.data <- rbind(plot.data, new.data)
}

if(color_metagenomics) {
  m <- ggplot(plot.data, aes(x=sample, y=indiv, color=metagenomics))
} else {
  m <- ggplot(plot.data, aes(x=sample, y=indiv))
}
  m <- m + geom_point() +
    theme_minimal() +
    xlab("") +
    ylab("") +
    theme(axis.text.x=element_blank()) +
    theme(axis.title.y=element_text(margin=margin(t=10, r=0, b=0, l=10)))
ggsave(paste0(plot_dir,"ABRP_overview.png"), scale=1.5, width=12, height=6, units="in")

}

# ====================================================================================================================
# sample frequency
# ====================================================================================================================

if(FALSE) {
# render histograms of sample frequency
histogram_indiv_samples(filtered)
histogram_sample_density(filtered, "weeks")
}

# ====================================================================================================================
# time courses for densely sampled individuals
# ====================================================================================================================

if(FALSE) {
short_list <- c("ACA", "DUI", "CAI", "COB", "DAS")
perform_mult_timecourse(filtered, short_list)
perform_mult_timecourse(data.piphillin, short_list)
}

# ====================================================================================================================
# covariance matrixies ordered on covariates
# ====================================================================================================================

# plot (subsampled) scrambled correlation matrices for a very coarse look into correlation on 
# the basis of various covariates
# batch stuff
visualize_groupwise_covariance(filtered, metadata, "plate", sample=50)                    # need to test these
visualize_groupwise_covariance(filtered, metadata, "flow_cell_lane", sample=1000)
visualize_groupwise_covariance(filtered, metadata, "library_pool", sample=200)
visualize_groupwise_covariance(filtered, metadata, "extract_dna_conc_ng", sample=100)
# biological stuff
visualize_groupwise_covariance(filtered, metadata, "sex", sample=1000)
visualize_groupwise_covariance(filtered, metadata, "season", sample=1000)
visualize_groupwise_covariance(filtered, metadata, "grp", sample=100)
visualize_groupwise_covariance(filtered, metadata, "matgrp", sample=100)
visualize_groupwise_covariance(filtered, metadata, "sname", sample=100)
visualize_groupwise_covariance(filtered, metadata, "age", sample=500)

# batch stuff
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "plate", sample=50)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "flow_cell_lane", sample=1000)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "library_pool", sample=200)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "extract_dna_conc_ng", sample=100)
# biological stuff
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "sex", sample=1000)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "season", sample=1000)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "grp", sample=100)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "matgrp", sample=100)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "sname", sample=100)
visualize_groupwise_covariance(data.piphillin, metadata.metagenomics, "age", sample=500)

# ====================================================================================================================
# ordination
# ====================================================================================================================

if(FALSE) {
# subset to a manageable sample size
md <- read_metadata(filtered)
sample_ids <- md$sample_id[sample(nsamples(filtered))[1:5000]]
baboon_counts <- subset_samples(filtered, sample_id %in% sample_ids)
# PCoA ordination using Bray-Curtis distance (phyloseq)
ord <- ordinate(baboon_counts, "PCoA", "bray")

# label SEASON
cat("Plotting seasonal labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="season")
ggsave(paste0(plot_dir,"season_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)

# label GROUP
# fix these categorical variables, so they don't accidentally plot as continuous values
cat("Plotting collection group-wise labeling on ordination...\n")
sample_data(baboon_counts)$grp <- as.factor(sample_data(baboon_counts)$grp)
p <- plot_ordination(baboon_counts, ord, type="sample", color="grp")
ggsave(paste0(plot_dir,"grp_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)

# label SEX
cat("Plotting sex-wise labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="sex")
ggsave(paste0(plot_dir,"sex_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)

# label MATERNAL GROUP
cat("Plotting maternal group-wise labeling on ordination...\n")
sample_data(baboon_counts)$matgrp <- as.factor(sample_data(baboon_counts)$matgrp)
p <- plot_ordination(baboon_counts, ord, type="sample", color="matgrp")
ggsave(paste0(plot_dir,"matgrp_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)

# label AGE
cat("Plotting maternal age-wise labeling on ordination...\n")
sample_data(baboon_counts)$age <- cut(sample_data(baboon_counts)$age, 
	breaks=c(0, 4.5, 19, Inf), labels=c("juvenile", "adult", "geriatric"))
p <- plot_ordination(baboon_counts, ord, type="sample", color="age")
ggsave(paste0(plot_dir,"age_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)

# resample where we want particularly samples witihn individuals
baboon_counts <- subset_samples(filtered, sname %in% best_sampled)
ord <- ordinate(baboon_counts, "PCoA", "bray")
# label INDIVIDUAL
cat("Plotting individual labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="sname")
ggsave(paste0(plot_dir,"sname_ordination_100.png"), plot=p, width=6, height=5, units="in", scale=1.5)
}
