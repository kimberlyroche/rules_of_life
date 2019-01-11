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
	#ilr_f1 <- apply_ilr(f1)
	#histogram_abundances(ilr_f1, "histogram_cutoff_3-20")
}

if(FALSE) {
	# filter to samples containing at least 90% counts >= 3 and get zero counts
	f2 <- filter_counts(glom_data, 3, 0.9)
	cat("Zero counts (90% filtering):",(1 - get_tiny_counts(f2, 1)),"\n")
	#ilr_f2 <- apply_ilr(f2)
	#histogram_abundances(ilr_f2, "histogram_cutoff_3-90")
}

if(TRUE) {
	baboons <- c("ACA", "DUI", "CAI", "COB", "DAS")
	for(b in baboons) {
		B_counts <- subset_samples(f2, sname==b)
		B_log_ratios <- apply_ilr(B_counts)
		plot_timecourse(B_counts, paste(b,"timecourse",sep="_"))
		# visualize correlation between times points (quick & dirty)
		# check (visually) preservation of collection_date order
		#check_sid_collection_order(B_counts)
		plot_corr_matrix(t(B_log_ratios), paste(b,"cor",sep="_"))
		# plot autocorrelation
		plot_autocorrelation(B_log_ratios,
				     read_metadata(B_counts),
				     lag.max=20,
				     filename=paste(b,"acf",sep="_"))
		plot_autocorrelation(B_log_ratios,
				     read_metadata(B_counts),
				     lag.max=20,
				     filename=paste(b,"acf","conservative",sep="_"),
				     norm_by_T=TRUE)
	}
}

if(FALSE) {
	cat("Writing histograms of sampling frequency...\n")
	#histogram_indiv_samples(f2)
	histogram_sample_density(f2, units="weeks")
}
