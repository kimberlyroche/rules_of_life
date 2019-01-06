source("include.R")

if(FALSE) {
	data <- read_data()
	cat("Zero counts (original):",(1 - get_tiny_counts(data, 1)),"\n")
	# remove technical replicates (for now)
	data <- subset_samples(data, sample_status==0)
	cat("Agglomerating data...\n")
	glom_data <- glom_counts(data)
	save(glom_data, file="glom_data_g.RData")
	#cat("Checking agglomeration...\n")
	#check_taxa_aggomeration(data)
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
	clr_f1 <- apply_clr(f1)
	histogram_abundances(clr_f1, "histogram_cutoff_3-20")
}

if(FALSE) {
	# filter to samples containing at least 90% counts >= 3 and get zero counts
	f2 <- filter_counts(glom_data, 3, 0.9)
	cat("Zero counts (90% filtering):",(1 - get_tiny_counts(f2, 1)),"\n")

	clr_f2 <- apply_clr(f2)
	histogram_abundances(clr_f2, "histogram_cutoff_3-90")
}

if(FALSE) {
	# filter to individuals "DUI" and "ACA"
	ACA_counts <- subset_samples(f2, sname=="ACA")
	ACA_log_ratios <- subset_samples(clr_f2, sname=="ACA")
	plot_timecourse(ACA_counts, "ACA_timecourse")
	# visualize correlation between times points (quick & dirty)
	ACA_data <- otu_table(ACA_log_ratios)@.Data
	# check (visually) preservation of collection_date order
	check_sid_collection_order(ACA_log_ratios)
	plot_corr_matrix(ACA_data, "ACA_cor")
	# plot autocorrelation
	plot_autocorrelation(ACA_log_ratios, lag.max=25, "ACA_acf")


	DUI_counts <- subset_samples(f2, sname=="DUI")
	DUI_log_ratios <- subset_samples(clr_f2, sname=="DUI")
	plot_timecourse(DUI_counts, "DUI_timecourse")
	# visualize correlation between times points (quick & dirty)
	DUI_data <- otu_table(DUI_log_ratios)@.Data
	# check (visually) preservation of collection_date order
	check_sid_collection_order(DUI_log_ratios)
	plot_corr_matrix(DUI_data, "DUI_cor")
	# plot autocorrelation
	plot_autocorrelation(DUI_log_ratios, lag.max=25, "DUI_acf")
}

if(TRUE) {
	cat("Writing histograms of sampling frequency...\n")
	histogram_indiv_samples(f2)
	histogram_sample_density(f2)
}
