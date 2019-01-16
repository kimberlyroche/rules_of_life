source("include.R")
library(driver)
#library(psych)
library(LaplacesDemon)

plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  ggsave(paste("plots/",filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# pull sample data for two individuals (DUI, ACA)

#load("glom_data_g.RData")
#f2 <- filter_counts(glom_data, 3, 0.9)

ilr_data <- NULL
week_kernel <- NULL
individual_kernel <- NULL
season_vector <- NULL

baboons <- c("ACA", "DUI", "CAI", "COB", "DAS")
for(b in 1:length(baboons)) {
  cat("Building kernels for baboon",baboons[b],"\n")
  baboon_counts <- subset_samples(f2, sname==baboons[b])
  baboon_metadata <- read_metadata(baboon_counts)
  baboon_log_ratios <- apply(otu_table(baboon_counts)@.Data+0.65, 1, ilr)

  nt <- ntaxa(baboon_counts)
  ns <- nsamples(baboon_counts)

  # build time kernel (distance per week)

  rho <- 1/2
  subset_week_kernel <- matrix(0, nrow=ns, ncol=ns)
  for(i in 1:ns) {
    for(j in 1:ns) {
      d1 <- as.Date(baboon_metadata$collection_date[i])
      d2 <- as.Date(baboon_metadata$collection_date[j])
      week_d <- abs(round((d2-d1)/7)) + 1
      # distance of 1 is < 7 days
      # distance of 2 is < 14 days, etc.
      if(week_d == 0) {
        week_d <- 1 # treat same day sampling as < 1 week
      }
      week_d <- week_d[[1]]
      subset_week_kernel[i,j] <- rho^week_d
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
  plot_cov(week_kernel, "date_correlation_test")

  indiv_cov <- 0.33
  subset_indiv_kernel <- matrix(0, nrow=ns, ncol=ns)
  subset_indiv_kernel[1:ns,1:ns] <- diag(ns)*(1-indiv_cov) + indiv_cov
  if(is.null(individual_kernel)) {
    individual_kernel <- subset_indiv_kernel
  } else {
    prev_r <- dim(individual_kernel)[1]
    prev_c <- dim(individual_kernel)[2]
    individual_kernel <- cbind(individual_kernel, matrix(0, nrow=dim(individual_kernel)[1], ncol=ns))
    individual_kernel <- rbind(individual_kernel, matrix(0, nrow=ns, ncol=dim(individual_kernel)[2]))
    new_r <- dim(individual_kernel)[1]
    new_c <- dim(individual_kernel)[2]
    individual_kernel[(prev_r+1):new_r,(prev_c+1):new_c] <- subset_indiv_kernel
  }
  plot_cov(individual_kernel, "indiv_correlation_test")

  subset_season_vector <- numeric(ns)
  for(i in 1:ns) {
    if(baboon_metadata$season[i] == "Wet") {
      subset_season_vector[i] <- 1
    } else {
      subset_season_vector[i] <- -1
    }
  }
  if(is.null(season_vector)) {
    season_vector <- subset_season_vector
  } else {
    season_vector <- c(season_vector, subset_season_vector)
  }
}

ns_all <- length(season_vector)
season_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
season_cov <- 0.3
for(i in 1:ns_all) {
  for(j in 1:ns_all) {
    if(i == j) {
      season_kernel[i,j] <- 1
    } else {
      season_kernel[i,j] <- season_vector[i]*season_vector[j]*season_cov
    }
  }
}

plot_cov(season_kernel, "season_correlation_test")
