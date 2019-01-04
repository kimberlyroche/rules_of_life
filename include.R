library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)
library(scales)
library(coda)

pdf(NULL)

# ====================================================================================================================
# ACCESSORS
# ====================================================================================================================

read_data <- function(write_sample=FALSE) {
	# rows are samples, columns are taxa
	cat("Reading data...\n")
	data <- readRDS("data/emp_baboon.RDS")
	if(write_sample) {
		write.table(otu_table(data)[1:10,1:10], file="otu_sample.txt", sep="\t")
		write.table(tax_table(data)[1:10,], file="tax_sample.txt", sep="\t")
	}
	return(data)
}

read_metadata <- function(data, write_sample=FALSE) {
	# metadata: rows are samples, "sample-specific variables" are columns
	cat("Reading metadata...\n")
	metadata <- phyloseq::sample_data(data)
	# sort ASC by collection date
	metadata <- metadata[order(metadata$sname,metadata$collection_date),]
	if(write_sample) {
		write.table(sample_data(data)[1:10,1:10], file="samp_sample.txt", sep="\t")
	}
	return(metadata)
}

# returns the (sampled) proportion at/above the given threshold
# so percent zeroes is 1 - get_tiny_counts(1)
get_tiny_counts <- function(data, threshold) {
	# just sample the samples
	use_ns <- 500
	if(nsamples(data) < 500) {
		use_ns <- nsamples(data)
	}
	sids <- sample(nsamples(data))[1:use_ns]
	return(sum(apply(otu_table(data)@.Data[sids,], c(1,2), function(x) { x >= threshold } ))/(ntaxa(data)*use_ns))
}

sid_to_collection_date <- function(data, sids) {
	metadata <- phyloseq::sample_data(data)
	metadata <- metadata[metadata$sample_id %in% sids,"collection_date"]
	return(unlist(metadata@.Data))
}

# check that order is preserved between sids >> collection dates
check_sid_collection_order <- function(data) {
	counts <- otu_table(data)@.Data
	sids <- dimnames(counts)[1][[1]]
	dates <- sid_to_collection_date(data, sids)
	print(order(dates))
}

check_taxa_aggomeration <- function(data) {
	# collapse on family produces 446 elements; check this is how many unique genus' there are
	tt <- apply(tax_table(data), 1, function(x) { paste(x[1:5], collapse="/") })
	tt_no <- unique(as.vector(tt))
	cat("Should == 446:",length(tt_no),"\n")
}

# ====================================================================================================================
# DATA TRANSFORMATION
# ====================================================================================================================

# uses phyloseq::filter_taxa to filter to taxa above a given count_threshold in at least 
# freq_threshold observations
filter_counts <- function(data, count_threshold, freq_threshold) {
	total_counts <- sum(otu_table(data)@.Data)
	filtered_data <- filter_taxa(data, function(x) sum(x >= count_threshold)/nsamples(data) >= freq_threshold, TRUE)
	total_counts_filtered <- sum(otu_table(filtered_data)@.Data)
	cat("Total counts:",total_counts,"filtered to:",total_counts_filtered,"\n")
	cat("Total taxa:",ntaxa(data),"filtered to:",ntaxa(filtered_data),"\n")
	return(filtered_data)
}

# NOTE: shouldn't agglomerate per individual but for testing, we'll try
glom_counts <- function(data, level="genus", NArm=FALSE) {
	cat("Agglomerating data to",level,"level...\n")
	glom_data <- tax_glom(data, taxrank=level, NArm=NArm)
	cat("Agglomerated data has dimension ntaxa=",ntaxa(glom_data),"x nsamples=",nsamples(glom_data),"\n")
	return(glom_data)
}

apply_proportion <- function(data) {
	return(transform_sample_counts(data, function(x) x / sum(x) ))
}

apply_clr <- function(data) {
	return(transform_sample_counts(data, function(x) { temp <- log(x+0.65); names(temp) <- names(x); temp; } ))
}

calc_autocorr <- function(data, lag.max=NULL) {
	# expects samples along rows, taxa along columns
	# this should be the correct average to take difference with respect to?
	# Prado & West pg. 9-10
	y.bar <- as.vector(apply(data, 2, mean)) # mean composition over time
	max.T <- dim(data)[1]
	if(is.null(lag.max)) {
		lag.max <- max.T-1
	}
	ac <- numeric(lag.max+1)
	gamma.0 <- -1
	for(lag in 0:lag.max) {
		tot <- 0
		for(t in 1:(max.T-lag)) {
			y.th <- as.vector(data[t+lag,])
			y.t <- as.vector(data[t,])
			tot <- tot + (t(y.th - y.bar)%*%(y.t - y.bar))[1]
		}
		gamma.h <- tot/max.T
		if(lag == 0) {
			gamma.0 <- gamma.h
		}
		ac[lag+1] <- gamma.h/gamma.0
	}
	return(ac)
}

get_num_indiv <- function(metadata) {
	individuals <- metadata$sname
	return(length(unique(individuals)))
}

# ====================================================================================================================
# VISUALIZATION
# ====================================================================================================================

# for all taxa, plots the percent at/above the given threshold, ordered descending
plot_percent_threshold <- function(data, threshold=3, save_filename) {
	cat("Plotting percent counts >=",threshold,"in all ASVs...\n")
	count_table <- otu_table(data)
	min_counts <- apply(count_table, 2, function(x) { sum(x >= threshold)/nsamples(data) })
	min_counts <- stack(sort(min_counts, decreasing=TRUE))
	p <- min_counts %>% ggplot(aes(x=seq(1,ntaxa(data)), y=values)) +
		geom_point(size=1) +
		theme_minimal() +
		xlab("ASV no.") +
		ylab(paste("Percent counts >= ",threshold,sep=""))
	ggsave(save_filename, plot=p, scale=1.5, width=5, height=3, units="in")
}

histogram_abundances <- function(data, filename="histogram") {
	df <- psmelt(data)
	p <- ggplot(df) +
		geom_histogram(aes(x=Abundance), binwidth=0.25) +
		theme_minimal() +
		xlab("log ratio abundance")
	ggsave(paste(filename,".png",sep=""), plot=p)
}

plot_corr_matrix <- function(data, filename) {
	cov_mat <- cov(t(data))
	cor_mat <- cov2cor(cov_mat)
	df <- gather_array(cor_mat)
	p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
		geom_tile() +
		xlab("samples") +
		ylab("samples") +
		theme_minimal() +
		guides(fill=guide_legend(title="correlation"))
	ggsave(paste(filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

plot_timecourse <- function(data, save_filename) {
	df <- psmelt(data)
	df2 <- bind_cols(list(OTU=df$OTU, Sample=df$Sample, Abundance=df$Abundance))

	# replace Sample ID's with their dates for readability
	metadata <- sample_data(data)
	for(i in 1:dim(df2)[1]) {
		df2$Sample[i] <- metadata[metadata$sample_id==df2$Sample[i],"collection_date"][[1]]
	}

	# make sure the dates are in order and fix the order by converting to factors
	df3 <- df2[order(df2$Sample),]
	df3$Sample <- factor(df3$Sample, levels=unique(df3$Sample))

	# replace ASV sequences with their (abbrev.) taxonomy for readability
	for(i in 1:dim(df3)[1]) {
		# show labels as order/family/genus
		# species is NA for all
		df3$OTU[i] <- paste(as.vector(tax_table(p)[df3$OTU[i],4:6]),collapse="/")
	}

	p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
		geom_bar(position="fill", stat="identity") +
		scale_y_continuous(labels = percent_format()) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		theme(legend.position="bottom") +
		theme(legend.text=element_text(size=8))
	ggsave(paste(save_filename,".png",sep=""), plot=p, scale=2, width=10, height=4, units="in")
}

histogram_indiv_samples <- function(data) {
	nt <- ntaxa(data)
	per_indiv <- psmelt(data) %>% group_by(sname) %>% summarise(n = n()/nt)
	p <- ggplot(per_indiv) +
		geom_histogram(aes(x=n), binwidth=5) +
		theme_minimal() +
		xlab("per-individual samples")
	ggsave("histogram_per_individual_samples.png", plot=p)
}

histogram_sample_density <- function(data) {
	per_indiv <- psmelt(data)
	snames <- unique(per_indiv$sname)
	differences <- c()
	min_diff <- Inf
	max_diff <- -Inf
	for(s in snames) {
		indiv_sname <- s
		df <- per_indiv[per_indiv$sname==indiv_sname,c("sample_id", "collection_date")]
		df <- unique(df)
		df <- df[order(df$collection_date),"collection_date"]
		if(length(df) > 1) {
			cat("Getting sample differences for individual",s,"...\n")
			indiv_differences <- as.numeric(length(df)-1)
			for(i in 1:(length(df)-1)) {
				d1 <- strptime(df[i], format="%Y-%m-%d")
				d2 <- strptime(df[i+1], format="%Y-%m-%d")
				diff <- d2-d1
				units(diff) <- "days"
				indiv_differences[i] <- as.numeric(diff)
				if(indiv_differences[i] > max_diff) {
					max_diff <- indiv_differences[i]
					cat("New MAX diff",max_diff,"from individual",s,"\n")
				}
				if(indiv_differences[i] < min_diff) {
					min_diff <- indiv_differences[i]
					cat("New MIN diff",min_diff,"from individual",s,"\n")
				}
			}
			differences <- c(differences, indiv_differences)
		}
	}
	diff_df <- as.data.frame(x=differences)
        p <- ggplot(diff_df) +
                geom_histogram(aes(x=differences), binwidth=14) +
                theme_minimal() +
                xlab("sample distance (days)")
        ggsave("histogram_sample_distance.png", plot=p)
	return(differences)
}

plot_autocorrelation <- function(data, lag.max, filename) {
	ac <- calc_autocorr(data, lag.max=lag.max)
	p <- ggplot(as.data.frame(cbind(x=seq(0,lag.max), y=ac)), aes(x=x, y=y)) +
		geom_point() +
		xlab("lag") +
		ylab("ACF") +
		theme_minimal()
	ggsave(paste(filename,".png",sep=""), plot=p, scale=1.5, width=3, height=3, units="in")
}
