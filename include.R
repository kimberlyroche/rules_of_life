library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)
library(scales)
library(coda)
library(Rcpp)

sourceCpp("fastCorr.cpp")

pdf(NULL)

best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")

# ====================================================================================================================
# ACCESSORS
# ====================================================================================================================

read_data <- function(write_sample=FALSE) {
  # rows are samples, columns are taxa
  data <- readRDS("data/emp_baboon.RDS")
  if(write_sample) {
    write.table(otu_table(data)[1:10,1:10], file="otu_sample.txt", sep="\t")
    write.table(tax_table(data)[1:10,], file="tax_sample.txt", sep="\t")
  }
  return(data)
}

read_metadata <- function(data, write_sample=FALSE) {
  # metadata: rows are samples, "sample-specific variables" are columns
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

perform_agglomeration <- function(level="genus") {
  data <- read_data()
  cat("Zero counts (original):",(1 - get_tiny_counts(data, 1)),"\n")
  # remove technical replicates (for now)
  data <- subset_samples(data, sample_status==0)
  cat("Agglomerating data...\n")
  glom_data <- glom_counts(data, level=level)
  save(glom_data, file=paste("glom_data_",level,".RData",sep=""))
  cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
  #cat("Checking agglomeration...\n")
  #check_taxa_aggomeration(data)
}

filter_data <- function(count_threshold=3, sample_threshold=0.9) {
  load("glom_data_genus.RData")
  cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
  filtered <- filter_counts(glom_data, count_threshold, sample_threshold)
  cat("Zero counts (filtered data):",(1 - get_tiny_counts(filtered, 1)),"\n")
  return(filtered)
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

# takes a phyloseq object, returns a data.frame!
apply_ilr <- function(data) {
  counts <- otu_table(data)@.Data
  return(apply(counts+0.65, 1, ilr))
}

# ====================================================================================================================
# VISUALIZATION -- HISTOGRAMS
# ====================================================================================================================

histogram_abundances <- function(data, filename="histogram") {
  df <- gather_array(data)
  p <- ggplot(df) +
    geom_histogram(aes(x=var), binwidth=0.25) +
    theme_minimal() +
    xlab("log ratio abundance")
  ggsave(paste(filename,".png",sep=""), plot=p)
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

histogram_sample_density <- function(data, units="days") {
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
        units(diff) <- units
        indiv_differences[i] <- as.numeric(diff)
        if(indiv_differences[i] > max_diff) {
          max_diff <- indiv_differences[i]
          cat("New MAX diff",max_diff,"from individual",s,"\n")
          print(d1)
          print(d2)
        }
        if(indiv_differences[i] < min_diff) {
          min_diff <- indiv_differences[i]
          cat("New MIN diff",min_diff,"from individual",s,"\n")
          print(d1)
          print(d2)
        }
      }
      differences <- c(differences, indiv_differences)
    }
  }
  diff_df <- as.data.frame(x=differences)
        p <- ggplot(diff_df) +
                geom_histogram(aes(x=differences), binwidth=4) +
                theme_minimal() +
                xlab(paste("sample distance (",units,")",sep=""))
        ggsave("histogram_sample_distance.png", plot=p)
}

# takes a list of correlations or covariances
histogram_corr <- function(data, filename, cov=FALSE) {
  df <- as.data.frame(data)
  colnames(df) <- c("corr")
  p <- ggplot(df) +
         geom_density(aes(x=corr)) +
         theme_minimal()
  if(cov) {
    p <- p + xlab("covariance")
  } else {
    p <- p + xlim(-1, 1) + xlab("correlation")
  }
  ggsave(paste(filename,".png",sep=""), plot=p)
}

# ====================================================================================================================
# VISUALIZATION -- COVARIANCE MATRICES/HEATMAPS
# ====================================================================================================================

plot_corr_matrix <- function(data, filename, cov=FALSE) {
  cor_mat <- fastCorr(t(data))
  #cov_mat <- cov(t(data))
  #cor_mat <- cov2cor(cov_mat)
  df <- gather_array(cor_mat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    scale_fill_gradientn(limits = c(-1,1),
    colours=c("navyblue", "darkmagenta", "darkorange1")) +
    xlab("samples") +
    ylab("samples") +
    theme_minimal()
  if(cov) {
    p <- p + guides(fill=guide_legend(title="covariance"))
  } else {
    p <- p + guides(fill=guide_legend(title="correlation"))
  }
  ggsave(paste(filename,".pdf",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

visualize_groupwise_covariance <- function(data, group, sample=1000000) {

  md <- read_metadata(data)

  if(group == "grp") {
    groups <- unique(md$grp); partition_obj <- md$grp; use_no <- sample
  } else if(group == "season") {
    groups <- unique(md$season); partition_obj <- md$season; use_no <- sample
  } else if(group == "matgrp") {
    groups <- unique(md$matgrp); partition_obj <- md$matgrp; use_no <- sample
  } else if(group == "sname") {
    #groups <- best_sampled; partition_obj <- md$sname; use_no <- sample
    groups <- unique(md$sname); partition_obj <- md$sname; use_no <- sample
  } else if(group == "sex") {
    groups <- unique(md$sex); partition_obj <- md$sex; use_no <- sample
  } else if(group == "age") {
    groups <- c(4.5, 19.0, 30.0); partition_obj <- md$age; use_no <- sample
  }

  partition_idx <- list()
  if(group != "age") {
    for(g in groups) {
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj==g, "sample_id"][[1]]
    }
  } else {
    for(g in 1:length(groups)) {
      gmin <- 0
      gmax <- groups[g]
      if(g > 1) {
        gmin <- groups[g-1]
      }
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj > gmin & partition_obj <= gmax,
                                                    "sample_id"][[1]]
    }
  }

  lr <- apply_ilr(data)
  lr <- t(apply(lr, 1, function(x) x - mean(x)))
  partition_samples <- list()
  for(i in 1:length(groups)) {
    partition_samples[[i]] <- lr[,partition_idx[[i]],drop=FALSE]
  }
  # test via unique(md[md$sample_id %in% partition_idx[[1]], "grp"]) etc.

  stacked_lr <- NULL
  samples <- list()
  for(i in 1:length(groups)) {
    total_samples <- dim(partition_samples[[i]])[2]
    use_ids <- sample(total_samples)
    use_upper <- use_no
    #if(total_samples < use_upper) {
    #  use_upper <- total_samples
    #}
    # omit groups with fewer than the requested number of samples
    if(total_samples >= use_upper) {
      cat(groups[i],"has",total_samples,"samples; using",use_upper,"\n")
      samples[[i]] <- partition_samples[[i]][,use_ids[1:use_upper]]
      if(is.null(stacked_lr)) {
        stacked_lr <- samples[[i]]
      } else {
        stacked_lr <- cbind(stacked_lr, samples[[i]])
      }
    }
  }
  cat("Sample matrix is",dim(stacked_lr)[1],"by",dim(stacked_lr)[2],"\n")

  plot_corr_matrix(t(stacked_lr), paste(group, "_cov_matrix", sep=""))
}

# ====================================================================================================================
# VISUALIZATION -- OTHER
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

plot_timecourse <- function(data, save_filename) {
  p <- apply_proportion(data)
  df <- psmelt(p)
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
    df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4:6]),collapse="/")
  }

  p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="bottom") +
    theme(legend.text=element_text(size=8))
  ggsave(paste(save_filename,".png",sep=""), plot=p, scale=2, width=10, height=4, units="in")
}

# baboons is a list of snames, e.g. c("ACA", "DUI", "CAI", "COB", "DAS")
perform_mult_timecourse <- function(data, baboons) {
  for(b in baboons) {
    b <<- b
    cat("Plotting timecourse for",b,"...\n")
    B_counts <- subset_samples(data, sname==b)
    B_log_ratios <- apply_ilr(B_counts)
    plot_timecourse(B_counts, paste(b[1],"timecourse",sep="_"))
    # check (visually) preservation of collection_date order
    #check_sid_collection_order(B_counts)
  }
}

plot_autocorrelation <- function(data, lag.max=26, filename="autocorrelation_all") {
  individuals <- unique(read_metadata(data)$sname)

  lag.sums <- numeric(lag.max)
  lag.measured <- numeric(lag.max)
  for(indiv in individuals) {
    # this weird syntactic hack seems to be necessary for subset_samples?
    # apparently the thing you're filtering against must be globally available
    indiv <<- indiv
    counts <- subset_samples(data, sname==indiv)
    log_ratios <- apply_ilr(counts)
    log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))
    cat(indiv[1],"has",dim(log_ratios)[2],"samples\n")
    md <- read_metadata(counts)
    # get distances between adjacent timepoints in weeks
    d1 <- as.Date(md$collection_date[1])
    week_d <- numeric(length(md$collection_date)-1)
    for(d in 2:length(md$collection_date)) {
      d2 <- as.Date(md$collection_date[d])
      week_d[d-1] <- round((d2-d1)/7) + 1
      d1 <- d2
    }
    for(lag in 1:lag.max) {
      idx <- which(week_d==lag)
      if(length(idx) >= 1) {
        tot <- 0
        for(t in idx) {
          y.t <- as.vector(log_ratios[,t])
          y.tt <- sqrt(y.t%*%y.t)
          y.h <- as.vector(log_ratios[,t+1])
          y.hh <- sqrt(y.h%*%y.h)
          tot <- tot + (y.t%*%y.h)/(y.tt*y.hh)
        }
        lag.measured[lag] <- lag.measured[lag] + length(idx)
        lag.sums[lag] <- lag.sums[lag] + tot
      }
    }
  }
  # print out autocorrelations
  for(lag in 1:lag.max) {
    cat("Lag",lag,"=",(lag.sums[lag]/lag.measured[lag]),"\n")
    lag.sums[lag] <- lag.sums[lag]/lag.measured[lag]
  }
  # plot
  p <- ggplot(as.data.frame(cbind(x=seq(1,lag.max), y=lag.sums)), aes(x=x, y=y)) +
    geom_line(linetype = "dashed") +
    scale_x_continuous(breaks=seq(0,lag.max,1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag (weeks)") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  ggsave("autocorrelation_all.png", plot=p, scale=2, width=4, height=3, units="in")
}
