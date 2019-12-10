source("include/R/general.R")

levels <- c("phylum", "family")
#levels <- c("genus")

if(FALSE) {
# HOW MANY INDIVIDUALS DOES THIS TAXON APPEAR IN >0 TIMES AT A >0 COUNT?
for(level in levels) {
	glom_data <- load_glommed_data(level)
	data <- subset_samples(glom_data, sname %in% sname_list) # 32 taxa x 9765 samples

	# HOW MANY THINGS ONLY APPEAR IN ONE INDIVIDUAL'S SAMPLES?
	min_count <- 1
	indiv_appearances <- numeric(ntaxa(data))
	for(i in 1:length(sname_list)) {
		indiv <- sname_list[i]
		counts <- otu_table(subset_samples(data, sname == indiv)) # samples x taxa
		indiv_appearances <- indiv_appearances +
			apply(counts, 2, function(x) { sum(sapply(x, function(y) y >= min_count)) > 1 } )
	}
	names(indiv_appearances) <- NULL
	cat("Level: ",toupper(level),"\n")
	cat("\t(At least",min_count,"in 1 sample...)\n")
	print(indiv_appearances)
	png(paste0("hist_",level,".png"))
	hist(indiv_appearances[indiv_appearances > 0],
		breaks=1:max(indiv_appearances))
	dev.off()
}
}

# HOW MANY TAXA ARE UNIQUE WITHIN THIS INDIVIDUAL?
for(level in levels) {
	glom_data <- load_glommed_data(level)
	data <- subset_samples(glom_data, sname %in% sname_list) # 32 taxa x 9765 samples

        threshold <- 5
	indiv_percents <- c()
	for(i in 1:length(sname_list)) {
		indiv <- sname_list[i]
                cat("Parsing",indiv,"\n")
		counts_in <- otu_table(subset_samples(data, sname == indiv)) # samples x taxa
		counts_out <- otu_table(subset_samples(data, sname != indiv)) # samples x taxa
		indiv_percents[i] <- sum(sapply(1:ntaxa(data), function(x) sum(counts_in[,x]) >= threshold & sum(counts_out[,x]) < threshold ))/ntaxa(data)
	}
	cat("Level: ",toupper(level),"\n")
        print(indiv_percents)
	png(paste0("hist_percent_",level,".png"))
	hist(indiv_percents,
		breaks=seq(0, 1, 0.01))
	dev.off()
}

