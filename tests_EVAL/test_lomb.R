library(lomb)
library(driver)
library(ggplot2)

level <- "family"
period <- 365

source("include.R")

glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.33, collapse_level=level, verbose=TRUE)

# repeat this for all individuals and see if any truly seasonal taxa jump out?

indiv_obj <- fitted_individuals(level="family")
#for(baboon in indiv_obj$individuals) {
for(baboon in c("WYN")) {
  cat("Individual",baboon,"\n")
  indiv_data <- subset_samples(data, sname==baboon)

  indiv_metadata <- read_metadata(indiv_data)
  baseline_date <- indiv_metadata$collection_date[1]
  observations <- sapply(indiv_metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1
  Y <- otu_table(indiv_data)@.Data
  clr.Y <- driver::clr(Y + 0.5) # samples x taxa

  for(i in 1:ncol(clr.Y)) {
    fit <- lsp(clr.Y[,i], times=observations, from=150, to=400, type="p", plot=TRUE)
    if(fit$p.value < 0.0001) {
      cat(i,":",fit$peak.at[1],",",fit$p.value,"\n")
      df <- data.frame(x=observations, y=clr.Y[,i])
      p <- ggplot(df, aes(x=x, y=y)) +
             geom_point() +
             geom_smooth(span=10)
      ggsave(p, file=paste0(baboon,"_",i,".png"), units="in", height=3, width=8, dpi=72, scale=1.5)
    }
  }
}
