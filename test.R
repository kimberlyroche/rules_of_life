library(stray)
library(driver)
library(ggplot2)

source("include/R/general.R")

pc <- 0.5
level <- "phylum"
other_idx <- 13

hosts <- over_40
hosts <- hosts[1:5]

phylum <- readRDS("data/filtered_phylum_5_20.rds")

df <- data.frame(time=c(), clrval=c(), host=c())

for(h in hosts) {
  cat("Host is:",h,"\n")
  phylum.subsampled <- subset_samples(phylum, sname==h)
  subsampled.md <- sample_data(phylum.subsampled)
  sids <- subsampled.md[order(subsampled.md$collection_date),]$sample_id
  dates <- subsampled.md[order(subsampled.md$collection_date),]$collection_date
  baseline <- dates[1]
  times <- sapply(dates, function(x) difftime(x, baseline, units="days")[[1]] )
  fitted.model <- readRDS(paste0(model_dir,level,"/",h,"_bassetfit.rds"))
  model.clr <- to_clr(fitted.model$fit)
  host.eta <- model.clr$Eta # dimensions: taxa x observations x posterior samples
                            # order of observations should be w/r/t time already
  posterior_samples <- gather_array(host.eta[other_idx,,], "other_value", "observation", "sample_no")
  posterior_samples$observation <- as.factor(posterior_samples$observation)
  levels(posterior_samples$observation) <- round(times)
  posterior_samples$observation <- as.numeric(posterior_samples$observation)

  post_quantiles <- posterior_samples %>%
    group_by(observation) %>%
    summarise(p2.5 = quantile(other_value, prob=0.025),
              p25 = quantile(other_value, prob=0.25),
              mean = mean(other_value),
              p75 = quantile(other_value, prob=0.75),
              p97.5 = quantile(other_value, prob=0.975)) %>%
    ungroup()

  p <- ggplot(post_quantiles, aes(x=observation, y=mean)) +
    geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    #geom_line(color="blue") +
    #geom_point(data=lr_tidy[lr_tidy$LR_coord==LR_coord,], aes(x=observation, y=LR_value), alpha=0.5) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45)) +
    ylab("LR coord")

  ggsave("test.png", plot=p, dpi=100, height=3, width=10, units="in")

  #df <- rbind(df, data.frame(time=times, clrval=phylum.clr[sids,other_idx], host=rep(h, length(times))))
}




# ===================================================================================================
# SAVED CODE
# ===================================================================================================

# these are pretty impossible to read; use estimates of eta from basset

# phylum.clr <- driver::clr(otu_table(phylum)@.Data + pc)

# for(h in hosts) {
#   cat("Host is:",h,"\n")
#   phylum.subsampled <- subset_samples(phylum, sname==h)
#   subsampled.md <- sample_data(phylum.subsampled)
#   sids <- subsampled.md[order(subsampled.md$collection_date),]$sample_id
#   dates <- subsampled.md[order(subsampled.md$collection_date),]$collection_date
#   baseline <- dates[1]
#   times <- sapply(dates, function(x) difftime(x, baseline, units="days")[[1]] )
#   df <- rbind(df, data.frame(time=times, clrval=phylum.clr[sids,other_idx], host=rep(h, length(times))))
# }

# p <- ggplot(df, aes(x=time, y=clrval)) +
#   geom_point() +
#   facet_wrap(vars(host), nrow=2, scales="free_x")
# ggsave("test.png", plot=p, dpi=100, height=4, width=10, units="in")
