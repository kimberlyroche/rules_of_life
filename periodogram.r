source("include.R")
library(lomb)

glom_data <- load_glommed_data(level="genus", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
# family-agglomerated has replicates; remove these
non_reps <- prune_samples(sample_data(filtered)$sample_status==0, filtered)
metadata <- read_metadata(non_reps)
cat("Filtered down to",nsamples(non_reps),"samples x",ntaxa(non_reps),"taxa\n")

sname <- "ACA"
sname <- "DUI"
indiv_samples <- prune_samples(sample_data(non_reps)$sname==sname, non_reps)
cat("Filtered down to",nsamples(indiv_samples),"samples x",ntaxa(indiv_samples),"taxa for individual",sname,"\n")

# expects samples as rows
counts <- otu_table(indiv_samples)@.Data
dates <- sample_data(indiv_samples)$collection_date
first_sample <- dates[1]
data.clr <- driver::clr(counts + 0.5)

df <- data.frame(date=c(0, sapply(dates[2:length(dates)], function(x) round(difftime(x, first_sample, units="days")))))
df <- cbind(df, data.clr)

for(i in 1:5) {
  fit <- lsp(df[,c(1,i)], from=7, to=730, type="period", alpha=0.05, plot=FALSE)
  plot(fit, type="l", main="Lomb-Scargle Periodogram", xlab=NULL, ylab="normalised power", level=FALSE, log="")
  plot(df[,1], df[,2], type="l", xlab="Days after first sample", ylab=paste0("CLR coord ",i))
}
