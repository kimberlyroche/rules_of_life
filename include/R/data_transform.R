base_path <- "/data/mukherjeelab/rulesoflife"

source(file.path(base_path,"include/R/general.R"))

# ====================================================================================================================
# DATA TRANSFORMATION, ETC.
# ====================================================================================================================

# recode date as numeric: number of days since earliest sample date
# these are frequently *fractional* because of leap years!!!
date_to_num <- function(date, baseline="2000-04-21") {
  as.numeric(difftime(date, as.Date(baseline), units="days"))
}

filter_data_omit <- function(count_threshold=3, sample_threshold=0.9, data=NULL) {
  if(is.null(data)) {
    load(file.path(data_dir,"glom_data_genus.RData"))
    #cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
    filtered <- filter_counts(glom_data, count_threshold, sample_threshold)
    #cat("Zero counts (filtered data):",(1 - get_tiny_counts(filtered, 1)),"\n")
  } else {
    filtered <- filter_counts(data, count_threshold, sample_threshold)
  }
  return(filtered)
}

# uses phyloseq::filter_taxa to filter to taxa above a given count_threshold in at least 
# freq_threshold observations
filter_counts <- function(data, count_threshold, freq_threshold) {
  total_counts <- sum(otu_table(data)@.Data)
  filtered_data <- filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= freq_threshold, TRUE)
  total_counts_filtered <- sum(otu_table(filtered_data)@.Data)
  #cat("Total counts:",total_counts,"filtered to:",total_counts_filtered,"\n")
  #cat("Total taxa:",ntaxa(data),"filtered to:",ntaxa(filtered_data),"\n")
  return(filtered_data)
}

# NOTE: shouldn't agglomerate per individual but for testing, we'll try
glom_counts <- function(data, level="genus", NArm=FALSE) {
  cat("Agglomerating data to",level,"level...\n")
  glom_data <- tax_glom(data, taxrank=level, NArm=NArm)
  cat("Agglomerated data has dimension ntaxa=",ntaxa(glom_data),"x nsamples=",phyloseq::nsamples(glom_data),"\n")
  return(glom_data)
}

apply_proportion <- function(data) {
  return(transform_sample_counts(data, function(x) x / sum(x) ))
}

# returns a data.frame
# if passing in a phyloseq object, expects taxa as columns
# if passing in metaganomics data, expects enzymes as rows
apply_ilr <- function(data, pseudocount=0.65) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data # samples (rows) x taxa (columns)
  } else {
    counts <- data # enzymes (rows) x samples (columns)
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  # driver log-ratio transforms expect samples as rows!
  # output: samples (rows) x D-1 taxa (columns)
  return(driver::ilr(counts))
}

# as above; returns D instead of D-1
apply_clr <- function(data, pseudocount=0.65) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data
  } else {
    counts <- data
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  return(driver::clr(counts))
}

# as above
# d is the index of the reference taxon
apply_alr <- function(data, pseudocount=0.65, d=NULL) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data
  } else {
    counts <- data
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  return(driver::alr(counts, d=d))
}

geom_mean <- function(x) {
  return(prod(x)**(1/length(x)))
}

aitchison_dist <- function(x.i, x.j) {
  sq.sum <- 0
  gm.i <- geom_mean(x.i)
  gm.j <- geom_mean(x.j)
  for(k in 1:length(x.i)) {
    sq.sum <- sq.sum + (log(x.i[k]/gm.i) - log(x.j[k]/gm.j))**2
  }
  return(sqrt(sq.sum))
}
