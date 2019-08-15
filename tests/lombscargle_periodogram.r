library(lomb)

source("../include/R/general.R")

level <- "genus"
data <- load_and_filter(level)
non_reps <- prune_samples(sample_data(data)$sample_status==0, data)
metadata <- read_metadata(non_reps)
cat("Filtered down to",nsamples(non_reps),"samples x",ntaxa(non_reps),"taxa\n")

plot_sample_periodogram <- TRUE
all_peaks <- numeric()
signif_peaks <- numeric()
for(sname in best_sampled) {
  indiv_samples <- prune_samples(sample_data(non_reps)$sname==sname, non_reps)
  cat("Filtered down to",nsamples(indiv_samples),"samples x",ntaxa(indiv_samples),"taxa for individual",sname,"\n")

  # expects samples as rows
  counts <- otu_table(indiv_samples)@.Data
  dates <- sample_data(indiv_samples)$collection_date
  first_sample <- dates[1]
  data.clr <- driver::clr(counts + 0.5)

  df <- data.frame(date=c(0, sapply(dates[2:length(dates)], function(x) round(difftime(x, first_sample, units="days")))))
  df <- cbind(df, data.clr)

  for(taxon in 1:ntaxa(indiv_samples)) {
    fit <- lsp(df[,c(1,taxon)], from=7, to=730, type="period", alpha=0.05, plot=FALSE)
    all_peaks[length(all_peaks)+1] <- fit$peak.at[1]
    if(fit$peak.at[2] <= 0.01) {
      signif_peaks[length(signif_peaks)+1] <- fit$peak.at[1]
      if(plot_sample_periodogram) {
        png(paste0(plot_dir,"periodogram_LS_sample_",sname,".png"), height=600, width=800)
        plot(fit, type="l", main=paste0("Lomb-Scargle Periodogram (",sname,")"), xlab=NULL, ylab="normalised power", level=FALSE, log="")
        dev.off()
        plot_sample_periodogram <- FALSE
      }
    }
  }
}

png(plot_dir,"periodogram_LS_allpeaks.png", height=600, width=800)
plot(density(all_peaks), main="Distribution of all peaks")
dev.off()
png(plot_dir,"periodogram_LS_signifpeaks.png", height=600, width=800)
plot(density(signif_peaks), main="Distribution of significant peaks (p <= 0.01)")
dev.off()
