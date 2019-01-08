source("include.R")

if(FALSE) {
	data <- read_data()
	cat("Zero counts (original):",(1 - get_tiny_counts(data, 1)),"\n")
	# remove technical replicates (for now)
	data <- subset_samples(data, sample_status==0)
	cat("Agglomerating data...\n")
	level <- "family"
	glom_data <- glom_counts(data, level=level)
	save(glom_data, file="glom_data_f.RData")
	cat("Checking agglomeration...\n")
	check_taxa_aggomeration(data)
}

if(FALSE) {
	load("glom_data_g.RData")
	# glom_data contains all samples agglomerated to genus level by phyloseq
	cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
}

if(FALSE) {
	# filter to samples containing at least 20% counts >= 3 and get zero counts
	f1 <- filter_counts(glom_data, 3, 0.2)
	cat("Zero counts (20% filtering):",(1 - get_tiny_counts(f1, 1)),"\n")

	# histogram log abundance for two zero-filtering levels
	ilr_f1 <- apply_ilr(f1)
	histogram_abundances(ilr_f1, "histogram_cutoff_3-20")
}

if(FALSE) {
	# filter to samples containing at least 90% counts >= 3 and get zero counts
	f2 <- filter_counts(glom_data, 3, 0.9)
	cat("Zero counts (90% filtering):",(1 - get_tiny_counts(f2, 1)),"\n")

	ilr_f2 <- apply_ilr(f2)
	histogram_abundances(ilr_f2, "histogram_cutoff_3-90")
}

if(TRUE) {
	# filter to individuals "DUI" and "ACA"
	ACA_counts <- subset_samples(f2, sname=="ACA")
	ACA_log_ratios <- apply_ilr(ACA_counts)
	plot_timecourse(ACA_counts, "ACA_timecourse")
	# visualize correlation between times points (quick & dirty)
	# check (visually) preservation of collection_date order
	#check_sid_collection_order(ACA_log_ratios)
	plot_corr_matrix(t(ACA_log_ratios), "ACA_cor")
	# plot autocorrelation
	plot_autocorrelation(ACA_log_ratios, lag.max=NULL, "ACA_acf")

	DUI_counts <- subset_samples(f2, sname=="DUI")
	DUI_log_ratios <- apply_ilr(DUI_counts)
	plot_timecourse(DUI_counts, "DUI_timecourse")
	#check_sid_collection_order(DUI_log_ratios)
	plot_corr_matrix(t(DUI_log_ratios), "DUI_cor")
	plot_autocorrelation(DUI_log_ratios, lag.max=NULL, "DUI_acf")
}

if(FALSE) {
	cat("Writing histograms of sampling frequency...\n")
	histogram_indiv_samples(f2)
	histogram_sample_density(f2)
}
