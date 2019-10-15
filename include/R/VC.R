source("include/R/general.R")
source("include/R/data_transform.R")
sourceCpp("include/cpp/dens_optim.cpp")

# ====================================================================================================================
# ESTIMATION OF VARIANCE COMPONENTS
# ====================================================================================================================

# plots crude variance components (per-individual) for diagnostic purposes
plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  ggsave(paste(plot_dir,filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# pass in zero-filtered data
estimate_variance_components <- function(data, metadata, optim_it=1, use_individuals=Inf, include_residual=False) {
  ilr_data <- NULL
  week_kernel <- NULL
  season_vector <- NULL
  age_vector <- NULL
  group_vector <- NULL
  indiv_vector <- NULL
  plate_vector <- NULL
  conc_vector <- NULL
  
  within_group_corr <- 0.1
  between_group_corr <- 0

  # build kernels for factors of interest

  baboons <- over_50
  if(use_individuals < Inf) {
    baboons <- baboons[sample(length(baboons))[1:use_individuals]] # subsample for debugging
  }
  for(b in 1:length(baboons)) {
    cat("Building kernel for",baboons[b],"\n")
    # get samples for this individual only
    if("phyloseq" %in% class(data)) {
      # phyloseq::subset_samples is causing all kinds of issues, so let's subset on sname manually!
      # basically you can't use subset_samples inside functions or loops
      # see: https://github.com/joey711/phyloseq/issues/752
      remove_idx <- as.character(get_variable(data, "sname")) == baboons[b]
      baboon_counts <- prune_samples(remove_idx, data)
      # subset metadata
      baboon_metadata <- read_metadata(baboon_counts)
    } else {
      baboon_counts <- subset_metagenomics_sname(data, baboons[b], metadata)
      # column names are sample IDs
      sample_ids <- colnames(baboon_counts)
      baboon_metadata <- metadata[metadata$sample_id %in% sample_ids,]
    }
    baboon_log_ratios <- t(apply_ilr(baboon_counts))
    # we'll ultimately pass ilr_data to optim ax taxa or enzymes (rows) x samples (columns)
    # so build it up in that orientation
    if(is.null(ilr_data)) {
      ilr_data <- baboon_log_ratios
    } else {
      # ilr_data is samples as rows, taxa or enzymes as columns
      ilr_data <- cbind(ilr_data, baboon_log_ratios)
    }

    if("phyloseq" %in% class(data)) {
      nt <- ntaxa(baboon_counts)
      ns <- phyloseq::nsamples(baboon_counts)
    } else {
      nt <- dim(baboon_log_ratios)[1] + 1 # no. taxa = no. enzymes; this is pre-LR transform
      ns <- dim(baboon_log_ratios)[2]
    }

    rho <- 1/2
    subset_week_kernel <- matrix(0, nrow=ns, ncol=ns)
    for(i in 1:ns) {
      for(j in 1:ns) {
        if(i == j) {
          subset_week_kernel[i,j] <- 1
        } else {
          d1 <- as.Date(baboon_metadata$collection_date[i])
          d2 <- as.Date(baboon_metadata$collection_date[j])
          # get number of 'weeks difference'
          week_d <- abs(round((d2-d1)/7))[[1]] + 1
          # distance of 1 is < 7 days
          # distance of 2 is < 14 days, etc.
          subset_week_kernel[i,j] <- rho^week_d + 0.1
        }
      }
    }

    if(is.null(week_kernel)) {
      week_kernel <- subset_week_kernel
    } else {
      prev_r <- dim(week_kernel)[1]
      prev_c <- dim(week_kernel)[2]
      week_kernel <- cbind(week_kernel, matrix(0, nrow=dim(week_kernel)[1], ncol=ns))
      week_kernel <- rbind(week_kernel, matrix(0, nrow=ns, ncol=dim(week_kernel)[2]))
      new_r <- dim(week_kernel)[1]
      new_c <- dim(week_kernel)[2]
      week_kernel[(prev_r+1):new_r,(prev_c+1):new_c] <- subset_week_kernel
    }

    subset_season_vector <- character(ns)
    subset_age_vector <- numeric(ns)
    subset_group_vector <- numeric(ns)
    subset_indiv_vector <- numeric(ns)
    subset_plate_vector <- numeric(ns)
    subset_conc_vector <- numeric(ns)
    for(i in 1:ns) {
      if(baboon_metadata$season[i] == "Wet") {
        subset_season_vector[i] <- "W"
      } else {
        subset_season_vector[i] <- "D"
      }
      if(baboon_metadata$age[i] < 4.5) {
        subset_age_vector[i] <- "J"
      } else if(baboon_metadata$age[i] > 19) {
        subset_age_vector[i] <- "G"
      } else {
        subset_age_vector[i] <- "A"
      }
      subset_group_vector[i] <- baboon_metadata$grp[i]
      subset_indiv_vector[i] <- baboon_metadata$sname[i]
      subset_plate_vector[i] <- baboon_metadata$plate[i]
      if(baboon_metadata$extract_dna_conc_ng[i] <= 2.5) {
        subset_conc_vector[i] <- 1
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 5.0) {
        subset_conc_vector[i] <- 2
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 7.5) {
        subset_conc_vector[i] <- 3
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 10.0) {
        subset_conc_vector[i] <- 4
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 15.0) {
        subset_conc_vector[i] <- 5
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 20.0) {
        subset_conc_vector[i] <- 6
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 25.0) {
        subset_conc_vector[i] <- 7
      } else {
        subset_conc_vector[i] <- 8
      }
    }
    if(is.null(season_vector)) {
      season_vector <- subset_season_vector
    } else {
      season_vector <- c(season_vector, subset_season_vector)
    }
    if(is.null(age_vector)) {
      age_vector <- subset_age_vector
    } else {
      age_vector <- c(age_vector, subset_age_vector)
    }
    if(is.null(group_vector)) {
      group_vector <- subset_group_vector
    } else {
      group_vector <- c(group_vector, subset_group_vector)
    }
    if(is.null(indiv_vector)) {
      indiv_vector <- subset_indiv_vector
    } else {
      indiv_vector <- c(indiv_vector, subset_indiv_vector)
    }
    if(is.null(plate_vector)) {
      plate_vector <- subset_plate_vector
    } else {
      plate_vector <- c(plate_vector, subset_plate_vector)
    }
    if(is.null(conc_vector)) {
      conc_vector <- subset_conc_vector
    } else {
      conc_vector <- c(conc_vector, subset_conc_vector)
    }
  }
  
  ns_all <- length(season_vector)

  season_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        season_kernel[i,j] <- 1
      } else if(season_vector[i] == "W" && season_vector[j] == "W") {
        season_kernel[i,j] <- within_group_corr
      } else if(season_vector[i] == "D" && season_vector[j] == "D") {
        season_kernel[i,j] <- within_group_corr
      } else {
        season_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  age_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        age_kernel[i,j] <- 1
      } else if(age_vector[i] == "J" && age_vector[j] == "J") {
        age_kernel[i,j] <- within_group_corr
      } else if(age_vector[i] == "A" && age_vector[j] == "A") {
        age_kernel[i,j] <- within_group_corr
      } else if(age_vector[i] == "G" && age_vector[j] == "G") {
        age_kernel[i,j] <- within_group_corr
      } else {
        age_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  group_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        group_kernel[i,j] <- 1
      } else if(group_vector[i] == group_vector[j]) {
        group_kernel[i,j] <- within_group_corr
      } else {
        group_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  indiv_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        indiv_kernel[i,j] <- 1
      } else if(indiv_vector[i] == indiv_vector[j]) {
        indiv_kernel[i,j] <- within_group_corr
      } else {
        indiv_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  plate_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        plate_kernel[i,j] <- 1
      } else if(plate_vector[i] == plate_vector[j]) {
        plate_kernel[i,j] <- within_group_corr
      } else {
        plate_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  conc_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        conc_kernel[i,j] <- 1
      } else if(conc_vector[i] == conc_vector[j]) {
        conc_kernel[i,j] <- within_group_corr
      } else {
        conc_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  # residual (white noise kernel)
  if(include_residual) {
    # not sure what I was going for here including group off-diagonal stuff
    #empty_kernel <- matrix(within_group_corr, ns_all, ns_all)
    #diag(empty_kernel) <- 1
    empty_kernel <- matrix(1, ns_all, ns_all)
  } else {
    empty_kernel <- matrix(0, ns_all, ns_all)
  }
  
  if(length(baboons) <= 0) {
    if("phyloseq" %in% class(data)) {
      tag <- ""
    } else {
      tag <- "metagenomics"
    }
    plot_cov(week_kernel, paste("date_sample_correlation", tag, sep=""))
    plot_cov(season_kernel, paste("season_sample_correlation", tag, sep=""))
    plot_cov(age_kernel, paste("age_sample_correlation", tag, sep=""))
    plot_cov(group_kernel, paste("group_sample_correlation", tag, sep=""))
    plot_cov(indiv_kernel, paste("indiv_sample_correlation", tag, sep=""))
    plot_cov(plate_kernel, paste("plate_sample_correlation", tag, sep=""))
    plot_cov(conc_kernel, paste("DNAconc_sample_correlation", tag, sep=""))
    if(include_residual) {
      plot_cov(empty_kernel, paste("residual_sample_correlation", tag, sep=""))
    }
  }

  # optimize the scale of each variance component
  logd_mat_t <- function(s) {
    return(logd_matrixt(s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8],
                        ilr_data, week_kernel, season_kernel, group_kernel, age_kernel,
                        indiv_kernel, plate_kernel, conc_kernel, empty_kernel))
  }
  
  components <- 8
  estimates <- matrix(0, components, optim_it)
  for(i in 1:optim_it) {
    cat("Optimization iteration",i,"\n")
    # about 3x faster if we have R optim call a compiled C++ function
    # could probably speed up further using DEoptim instead of optim here
    res <- optim(par=runif(components), fn=logd_mat_t, method="L-BFGS-B",
                 lower=rep(0.00001,components), upper=rep(10,components))
    estimates[,i] <- res$par
  }
  
  kernels <- c("weekly","seasonal","group","age","individual","plate","DNAconc","residual")
  if(optim_it > 1) {
    if(include_residual) {
      total_var <- colSums(estimates)
    } else {
      total_var <- colSums(estimates[1:7,])
    }
  } else {
    if(include_residual) {
      total_var <- c(sum(estimates))
    } else {
      total_var <- c(sum(estimates[1:7,]))
    }
  }
  for(i in 1:optim_it) {
    rank <- order(estimates[,i], decreasing=TRUE)
    cat("Iteration",i,":")
    for(j in 1:components) {
      if(rank[j] != components || include_residual) {
        cat(" ",kernels[rank[j]]," (",round(estimates[rank[j],i]/total_var[i],digits=3),")",sep="")
      }
    }
    cat("\n")
  }
}
