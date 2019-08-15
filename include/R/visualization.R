source("include/R/general.R")

# ====================================================================================================================
# VISUALIZATION -- HISTOGRAMS
# ====================================================================================================================

histogram_abundances <- function(data, filename="histogram") {
  df <- gather_array(data)
  p <- ggplot(df) +
    geom_histogram(aes(x=var), binwidth=0.25) +
    theme_minimal() +
    xlab("log ratio abundance")
  ggsave(paste0(plot_dir,filename,".png",sep=""), plot=p)
}

histogram_indiv_samples <- function(data) {
  nt <- ntaxa(data)
  per_indiv <- psmelt(data) %>% group_by(sname) %>% summarise(n = n()/nt)
  p <- ggplot(per_indiv) +
    geom_histogram(aes(x=n), binwidth=5) +
    theme_minimal() +
    xlab("per-individual samples")
  ggsave(paste0(plot_dir,"histogram_per_individual_samples.png"), plot=p)
}

histogram_sample_density <- function(data, units="weeks") {
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
  if(units == "weeks") {
    p <- p + xlim(0, 100)
  }
  ggsave(paste0(plot_dir,"histogram_sample_distance.png"), plot=p)
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
  ggsave(paste0(plot_dir,filename,".png"), plot=p)
}

# ====================================================================================================================
# VISUALIZATION -- COVARIANCE MATRICES/HEATMAPS
# ====================================================================================================================

plot_corr_matrix <- function(data, filename, cov=FALSE) {
  cor_mat <- fast_corr(t(data))
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
  ggsave(paste0(plot_dir,filename,".png"), plot=p, scale=1.5, width=4, height=3, units="in")
}

visualize_groupwise_covariance <- function(data, md, group, sample=1000000) {

  if(group == "grp") {
    groups <- unique(md$grp); partition_obj <- md$grp; use_no <- sample
  } else if(group == "season") {
    groups <- unique(md$season); partition_obj <- md$season; use_no <- sample
  } else if(group == "matgrp") {
    groups <- unique(md$matgrp); partition_obj <- md$matgrp; use_no <- sample
  } else if(group == "sname") {
    # subsample the individuals
    #groups <- best_sampled; partition_obj <- md$sname; use_no <- sample
    groups <- unique(md$sname); partition_obj <- md$sname; use_no <- sample
  } else if(group == "sex") {
    groups <- unique(md$sex); partition_obj <- md$sex; use_no <- sample
  } else if(group == "age") {
    groups <- c(4.5, 19.0, 30.0); partition_obj <- md$age; use_no <- sample
  } else if(group == "plate") {
    # subsample the plates
    groups <- sample(unique(md$plate))[1:20]; partition_obj <- md$plate; use_no <- sample
  } else if(group == "flow_cell_lane") {
    groups <- unique(md$flow_cell_lane); partition_obj <- md$flow_cell_lane; use_no <- sample
  } else if(group == "library_pool") {
    groups <- unique(md$library_pool); partition_obj <- md$library_pool; use_no <- sample
  } else if(group == "extract_dna_conc_ng") {
    groups <- c(2.5, 5, 7.5, 10, 15, 20, 25, 100); partition_obj <- md$extract_dna_conc_ng; use_no <- sample
  }
  
  # list corresponding to elements of groups with sample IDs in each group's bucket
  partition_idx <- list()
  if(group == "age" || group == "extract_dna_conc_ng") {
    for(g in 1:length(groups)) {
      gmin <- 0
      gmax <- groups[g]
      if(g > 1) {
        gmin <- groups[g-1]
      }
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj > gmin & partition_obj <= gmax,
                                                    "sample_id"][[1]]
    }
  } else {
    for(g in groups) {
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj==g, "sample_id"][[1]]
    }
  }

  lr <- apply_ilr(data) # after this data is samples x taxa or enzymes
  lr <- scale(lr, center=T, scale=F)
  partition_samples <- list()
  for(i in 1:length(groups)) {
    if("phyloseq" %in% class(data)) {
      partition_samples[[i]] <- lr[partition_idx[[i]],,drop=FALSE]
    } else {
      partition_samples[[i]] <- lr[rownames(lr) %in% partition_idx[[i]],,drop=FALSE]
    }
  }
  # test clean partition via unique(md[md$sample_id %in% partition_idx[[1]], "grp"]) etc.
  
  stacked_lr <- NULL
  samples <- list()
  for(i in 1:length(groups)) {
    total_samples <- dim(partition_samples[[i]])[1]
    use_ids <- sample(total_samples)
    use_upper <- use_no
    #if(total_samples < use_upper) {
    #  use_upper <- total_samples
    #}
    # omit groups with fewer than the requested number of samples
    if(total_samples >= use_upper) {
      cat(groups[i],"has",total_samples,"samples; using",use_upper,"\n")
      samples[[i]] <- partition_samples[[i]][use_ids[1:use_upper],]
      if(is.null(stacked_lr)) {
        stacked_lr <- samples[[i]]
      } else {
        #stacked_lr <- cbind(stacked_lr, samples[[i]])
        stacked_lr <- rbind(stacked_lr, samples[[i]])
      }
    }
  }
  cat("Sample matrix is",dim(stacked_lr)[1],"by",dim(stacked_lr)[2],"\n")
  if("phyloseq" %in% class(data)) {
    save_filename <- paste("plots/", group, "_cov_matrix", sep="")
  } else {
    save_filename <- paste("plots/", group, "_cov_matrix_metagenomics", sep="")
  }
  save_filename <- 
  plot_corr_matrix(stacked_lr, save_filename)
}

# ====================================================================================================================
# VISUALIZATION -- OTHER
# ====================================================================================================================

# for all taxa, plots the percent at/above the given threshold, ordered descending
plot_percent_threshold <- function(data, threshold=3, save_filename) {
  cat("Plotting percent counts >=",threshold,"in all ASVs...\n")
  count_table <- otu_table(data)
  min_counts <- apply(count_table, 2, function(x) { sum(x >= threshold)/phyloseq::nsamples(data) })
  min_counts <- stack(sort(min_counts, decreasing=TRUE))
  p <- min_counts %>% ggplot(aes(x=seq(1,ntaxa(data)), y=values)) +
    geom_point(size=1) +
    theme_minimal() +
    xlab("ASV no.") +
    ylab(paste("Percent counts >= ",threshold,sep=""))
  ggsave(paste0(plot_dir,save_filename,".png"), plot=p, scale=1.5, width=5, height=3, units="in")
}

# plot proportional change over time
# note: palette size would probably have to increase to accomodate legend_level finer than family!
plot_timecourse_phyloseq <- function(data, save_filename, gapped=FALSE,
                                     legend=TRUE, legend_level="family", selected_samples=NULL) {
  n <- phyloseq::nsamples(data)
  p <- apply_proportion(data)
  df <- psmelt(p)
  df2 <- bind_cols(list(OTU=df$OTU, Sample=df$Sample, Abundance=df$Abundance, SID=df$sid))
  if(!is.null(selected_samples)) {
    df2$alpha <- 0.25
    df2[df2$SID %in% selected_samples,]$alpha <- 1
  }

  # for gapped plots, this is the OTU ID placeholder we'll use for dummy data
  na.string <- ".N/A"
  
  # replace Sample ID's with their dates for readability
  metadata <- sample_data(data)
  for(i in 1:dim(df2)[1]) {
    df2$Sample[i] <- paste(metadata[metadata$sample_id==df2$Sample[i],"collection_date"][[1]],
                           metadata[metadata$sample_id==df2$Sample[i],"sid"][[1]])
  }

  # make sure the dates are in order and fix the order by converting to factors
  df3 <- df2[order(df2$Sample),]
  if(gapped) {
    # insert empty samples were gaps of 2 weeks or more exist
    gap.days <- 13
    dates_present <- unique(df3$Sample)
    for(d in 1:(length(dates_present)-1)) {
      diff <- as.Date(dates_present[d+1]) - as.Date(dates_present[d])
      next.date <- as.Date(dates_present[d])
      attr(diff, "units") <- "days"
      while(diff > gap.days) {
        next.date <- next.date + gap.days
        df3 <- rbind(df3, list(OTU=na.string, Sample=as.character(next.date), Abundance=1))
        diff <- as.Date(dates_present[d+1]) - next.date
      }
    }
  }
  df3 <- df3[order(df3$Sample),]

  df3$Sample <- factor(df3$Sample, levels=unique(df3$Sample))
  level_idx <- 7
  if(legend_level == "kingdom") { level_idx = 1; }
  if(legend_level == "phylum") { level_idx = 2; }
  if(legend_level == "class") { level_idx = 3; }
  if(legend_level == "order") { level_idx = 4; }
  if(legend_level == "family") { level_idx = 5; }
  if(legend_level == "genus") { level_idx = 6; }
  # replace ASV sequences with their (abbrev.) taxonomy for readability
  for(i in 1:dim(df3)[1]) {
    # show labels as order/family/genus
    # species is NA for all
    if(df3$OTU[i] != na.string) {
      df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],level_idx]),collapse="/")
    }
  }

  categories <- unique(df3$OTU)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df3$OTU)))
  if(gapped) {
    coul[1] <- "#DDDDDD"
    img_width <- 15
  } else {
    img_width <- 10
  }
  if(!is.null(selected_samples)) {
    p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU, alpha=alpha)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8)) +
      theme(legend.position="none")
  } else {
    p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8))
    if(legend) {
      p <- p + theme(legend.position="bottom")
    } else {
      p <- p + theme(legend.position="none")
    }
  }
  if(n < 20) {
    # these are likely replicates
    img_width <- 5
  }
  ggsave(paste0(save_filename,".png"), plot=p, scale=1.5, dpi=100, width=img_width, height=3, units="in")
}

# this is just a wrapper to call plot_timecourse_phyloseq on a representative sample of individuals
# baboons is a list of snames, e.g. c("ACA", "DUI", "CAI", "COB", "DAS")
perform_mult_timecourse <- function(data, baboons, gapped=FALSE, legend=FALSE) {
  for(b in baboons) {
    b <<- b
    cat("Plotting timecourse for",b,"...\n")
    B_counts <- subset_samples(data, sname==b)
    B_log_ratios <- apply_ilr(B_counts)
    if(gapped) {
      save_filename <- paste("plots/",b[1],"_gapped_timecourse",sep="")
    } else {
      save_filename <- paste("plots/",b[1],"_timecourse",sep="")
    }
    plot_timecourse_phyloseq(B_counts, save_filename=save_filename, gapped=gapped, legend=legend)
  }
}

calc_autocorrelation <- function(data,
                                 metadata,
                                 resample=FALSE,
                                 lag.max=26,
                                 date_diff_units="weeks",
                                 use_lr="clr",
                                 alr_ref=NULL,
                                 resample_rate=0.2) {
  if(date_diff_units != "weeks" && date_diff_units != "months" && date_diff_units != "seasons") {
    date_diff_units <- "weeks"
  }

  individuals <- unique(metadata$sname)
  #individuals <- individuals[1:20] # for testing

  # only if using 
  season_boundaries <- c(200005, 200010, 200105, 200110, 200205, 200210, 200305, 200310,
                       200405, 200410, 200505, 200510, 200605, 200610, 200705, 200710,
                       200805, 200810, 200905, 200910, 201005, 201010, 201105, 201110,
                       201205, 201210, 201305, 201310)

  rounds <- 1
  if(resample) {
    rounds <- 100
    #rounds <- 10 # for testing
  }
  lags <- matrix(0, nrow=lag.max, ncol=rounds)

  # TBD: try various CoDa proportionality measures
  if(use_lr == "clr") {
    log_ratios <- apply_clr(data)
  } else if(use_lr == "alr") {
    log_ratios <- apply_alr(data, d=alr_ref)
  } else {
    log_ratios <- apply_ilr(data)
  }
  log_ratios <- scale(log_ratios, center=TRUE, scale=FALSE)

  for(r in 1:rounds) {
    if(resample) {
      cat("Resampling iteration",r,"\n")
    }
    lag.sums <- numeric(lag.max)
    lag.measured <- numeric(lag.max)
    for(indiv in individuals) {
      cat("Individual",indiv,"...\n")
      # pick an individual
      # this weird syntactic hack seems to be necessary for subset_samples?
      # apparently the thing you're filtering against must be globally available
      indiv <<- indiv
      # pull sample IDs associated with this individual
      sample_info <- metadata[metadata$sname == indiv, c("sample_id", "collection_date")]
      # appear already ordered but let's be paranoid
      sample_info <- sample_info[order(sample_info$collection_date),]
      # we need individuals with at least 2 samples
      if(length(intersect(rownames(log_ratios), sample_info$sample_id)) > 1) {
        sample_lr <- log_ratios[rownames(log_ratios) %in% sample_info$sample_id,,drop=F]
        indiv_sample_no <- dim(sample_lr)[1]
        # randomly censor some of the samples
        do_sample <- rep(1, indiv_sample_no)
        if(resample) {
          do_sample <- rbinom(indiv_sample_no, 1, resample_rate)
        }
        # replace these sample IDs with collection dates, that's what we really want
        # order should be maintained here
        rownames(sample_lr) <- sample_info[sample_info$sample_id %in% rownames(sample_lr), "collection_date"]$collection_date
  
        # get distances between adjacent timepoints in units of {date_diff_units}
        d1 <- as.Date(sample_info$collection_date[1])
        time_diff <- matrix(0, nrow=indiv_sample_no, ncol=indiv_sample_no)
        if(indiv_sample_no > 1) {
          for(d1_idx in 1:(indiv_sample_no-1)) {
            for(d2_idx in (d1_idx+1):indiv_sample_no) {
              d1 <- as.Date(sample_info$collection_date[d1_idx])
              d2 <- as.Date(sample_info$collection_date[d2_idx])
              time_diff[d1_idx,d2_idx] <- as.numeric(difftime(d2, d1, units="weeks"))
              if(date_diff_units == "weeks") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx])
              } else if(date_diff_units == "months") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx]/4)
              } else if(date_diff_units == "seasons") {
                d1_ym <- as.numeric(format(d1, "%Y%m"))
                d2_ym <- as.numeric(format(d2, "%Y%m"))
                # distance is number of season boundaries between the samples, indicated in the vector above
                time_diff[d1_idx,d2_idx] <- sum(season_boundaries < d2_ym) - sum(season_boundaries < d1_ym) + 1
              }
            }
          }
        }
        for(lag in 1:lag.max) {
          # find pairs with this lag
          idx <- which(time_diff==lag, arr.ind=T)
          if(dim(idx)[1] >= 1) {
            for(t in 1:dim(idx)[1]) {
              d1_idx <- idx[t,1]
              d2_idx <- idx[t,2]
              if(do_sample[d1_idx] && do_sample[d2_idx]) {
                # given centering, cosine angle == Pearson's correlation
                lag.sums[lag] <- lag.sums[lag] + cor(sample_lr[d1_idx,], sample_lr[d2_idx,])
                lag.measured[lag] <- lag.measured[lag] + 1
              }
            }
          }
          if(lag.sums[lag] > 0 & lag.measured[lag] > 0) {
            lags[lag,r] <- lag.sums[lag]/lag.measured[lag]
          }
        }
      }
    }
    if(rounds == 1) {
      # print out autocorrelations
      for(lag in 1:lag.max) {
        cat("Lag",lag,"=",lags[lag,r],"\n")
      }
      cat("\n")
    }
  }

  if(resample) {
    # get 90% confidence intervals, summarize
    lower_CI <- as.numeric(lag.max)
    upper_CI <- as.numeric(lag.max)
    avg_ac <- as.numeric(lag.max)
    for(lag in 1:lag.max) {
      ac_measured <- sort(lags[lag,])
      lower_CI[lag] <- ac_measured[round(rounds*0.1)]
      upper_CI[lag] <- ac_measured[round(rounds*0.9)]
      avg_ac[lag] <- mean(ac_measured)
    }
    return(data.frame(cbind(lag=seq(1,lag.max), ac=avg_ac, lower=lower_CI, upper=upper_CI)))
  } else {
    return(as.data.frame(cbind(lag=seq(1,lag.max), ac=as.vector(lags))))
  }
}

plot_mean_autocorrelation <- function(lags, filename="autocorrelation", width=6, height=3) {
  lag.max <- max(lags$lag)
  p <- ggplot(lags, aes(x=lag, y=ac)) +
    geom_line(linetype = "dashed") +
    scale_x_continuous(breaks=seq(0,lag.max,1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  ggsave(paste0(plot_dir,filename,".png"), plot=p, dpi=100, scale=1.5, width=width, height=height, units="in")
}

plot_bounded_autocorrelation <- function(lags, filename="autocorrelation", width=6, height=3) {
  p <- ggplot(lags, aes(x=lag, y=ac)) +
    geom_line(linetype="solid") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    scale_x_continuous(breaks=seq(0,dim(lags)[1],1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  ggsave(paste0(plot_dir,filename,".png"), plot=p, dpi=100, scale=1.5, width=width, height=height, units="in")
}
