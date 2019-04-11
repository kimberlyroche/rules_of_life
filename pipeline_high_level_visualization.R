source("include.R")

# ====================================================================================================================
# plot all samples (time on x-axis, individual on y-axis)
# ====================================================================================================================

# recode date as numeric: number of days since earliest sample date
date_to_num <- function(date) {
  as.numeric(difftime(date, as.Date("2000-04-21"), units="days"))
}

# read raw data
data <- read_data()
read_metadata(data)

cat("Earliest sample:", min(metadata$collection_date), "\n")
cat("Latest sample:", max(metadata$collection_date), "\n")

min_samples <- 10 # drop individuals with fewer than this many samples
individuals <- unique(metadata$sname)
samples <- list()
for(i in 1:length(individuals)) {
  indiv_samples <- subset_samples(metadata, (sname %in% c(individuals[i])))$collection_date
  indiv_samples <- as.vector(sapply(indiv_samples, date_to_num))
  if(length(indiv_samples) >= min_samples) {
    samples[[individuals[i]]] <- indiv_samples
  }
}

plot.data <- data.frame(indiv=NULL, sample=NULL)
for(ind in names(samples)) {
  new.data <- data.frame(indiv=rep(ind, length(samples[[ind]])), sample=samples[[ind]])
  plot.data <- rbind(plot.data, new.data)
}

m <- ggplot(plot.data, aes(x=sample, y=indiv)) +
  geom_point() +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(margin=margin(t=10, r=0, b=0, l=10)))
ggsave("plots/ABRP_overview.png", scale=2, width=4, height=8, units="in")

# ====================================================================================================================
# sample frequency
# ====================================================================================================================

glom_data <- load_glommed_data(level="species", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)

# render histograms of sample frequency
histogram_indiv_samples(filtered)
histogram_sample_density(filtered, "weeks")

# ====================================================================================================================
# time courses for densely sampled individuals
# ====================================================================================================================

perform_mult_timecourse(filtered, c("ACA", "DUI", "CAI", "COB", "DAS"))

# ====================================================================================================================
# covariance matrixies ordered on covariates
# ====================================================================================================================

# plot (subsampled) scrambled correlation matrices for a very coarse look into correlation on 
# the basis of various covariates
# batch stuff
visualize_groupwise_covariance(filtered, "plate", sample=50)
visualize_groupwise_covariance(filtered, "flow_cell_lane", sample=1000)
visualize_groupwise_covariance(filtered, "library_pool", sample=200)
visualize_groupwise_covariance(filtered, "extract_dna_conc_ng", sample=100)
# biological stuff
visualize_groupwise_covariance(filtered, "sex", sample=1000)
visualize_groupwise_covariance(filtered, "season", sample=1000)
visualize_groupwise_covariance(filtered, "grp", sample=100)
visualize_groupwise_covariance(filtered, "matgrp", sample=100)
visualize_groupwise_covariance(filtered, "sname", sample=100)
visualize_groupwise_covariance(filtered, "age", sample=500)

# ====================================================================================================================
# ordination
# ====================================================================================================================

# subset to a manageable sample size
md <- read_metadata(filtered)
sample_ids <- md$sample_id[sample(nsamples(filtered))[1:5000]]
baboon_counts <- subset_samples(filtered, sample_id %in% sample_ids)
# PCoA ordination using Bray-Curtis distance (phyloseq)
ord <- ordinate(baboon_counts, "PCoA", "bray")

# label SEASON
cat("Plotting seasonal labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="season")
ggsave("plots/season_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

# label GROUP
# fix these categorical variables, so they don't accidentally plot as continuous values
cat("Plotting collection group-wise labeling on ordination...\n")
sample_data(baboon_counts)$grp <- as.factor(sample_data(baboon_counts)$grp)
p <- plot_ordination(baboon_counts, ord, type="sample", color="grp")
ggsave("plots/grp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

# label SEX
cat("Plotting sex-wise labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="sex")
ggsave("plots/sex_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

# label MATERNAL GROUP
cat("Plotting maternal group-wise labeling on ordination...\n")
sample_data(baboon_counts)$matgrp <- as.factor(sample_data(baboon_counts)$matgrp)
p <- plot_ordination(baboon_counts, ord, type="sample", color="matgrp")
ggsave("plots/matgrp_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

# label AGE
cat("Plotting maternal age-wise labeling on ordination...\n")
sample_data(baboon_counts)$age <- cut(sample_data(baboon_counts)$age, 
	breaks=c(0, 4.5, 19, Inf), labels=c("juvenile", "adult", "geriatric"))
p <- plot_ordination(baboon_counts, ord, type="sample", color="age")
ggsave("plots/age_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)

# resample where we want particularly samples witihn individuals
baboon_counts <- subset_samples(filtered, sname %in% best_sampled)
ord <- ordinate(baboon_counts, "PCoA", "bray")
# label INDIVIDUAL
cat("Plotting individual labeling on ordination...\n")
p <- plot_ordination(baboon_counts, ord, type="sample", color="sname")
ggsave("plots/sname_ordination_100.png", plot=p, width=6, height=5, units="in", scale=1.5)
