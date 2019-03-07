library(driver)
library(phyloseq)
library(ggplot2)

source("include.R")

translate_month <- function(x) {
  if(x == 2) { return("Feb") }
  #if(x == 3) { return("Mar") }
  if(x == 4) { return("Apr") }
  if(x == 5) { return("May") }
  if(x == 6) { return("Jun") }
  if(x == 7) { return("Jul") }
}

load("glom_data_genus_reps.RData")
md <- sample_data(glom_data)

# get correlation for replicates

# this should be agglomerated data
filtered <- filter_data(count_threshold=3, sample_threshold=0.2, data=glom_data)
counts <- otu_table(filtered)@.Data
log_ratios <- apply(counts + 0.65, 1, alr)
log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))

cat("Zero-filtered taxa to:",ntaxa(filtered),"\n")
replicate_samples <- md[md$sample_status==2,]
unique_replicates <- unique(replicate_samples[,c("sname","collection_date")])
corr_list <- c()
min_corr <- c(Inf, "", "")
max_corr <- c(-Inf, "", "")
for(i in 1:dim(unique_replicates)[1]) {
  sname <- unique_replicates[i]$sname
  date <- unique_replicates[i]$collection_date
  idx <- as.character(get_variable(replicate_samples, "sname")) == sname
  samples <- prune_samples(idx, replicate_samples)
  idx <- as.character(get_variable(samples, "collection_date")) == date
  samples <- prune_samples(idx, samples)
  # mean-center
  sample_lr <- log_ratios[,colnames(log_ratios) %in% samples$sample_id]
  # calculate correlation
  # the plots of log counts against each other look good
  # why is the "correlation" so low or even negative?
  total_corr <- 0
  pairs <- 0
  for(j in 1:(dim(sample_lr)[2]-1)) {
    for(k in 2:dim(sample_lr)[2]) {
      if(j != k) {
        y.t <- as.vector(sample_lr[,j])
        y.tt <- sqrt(y.t%*%y.t)
        y.h <- as.vector(sample_lr[,k])
        y.hh <- sqrt(y.h%*%y.h)
        #cat("Corr:",((y.t%*%y.h)/(y.tt*y.hh)),"\n")
        total_corr <- total_corr + (y.t%*%y.h)/(y.tt*y.hh)
        pairs <- pairs + 1
      }
    }
  }
  corr_list[i] <- total_corr/pairs
  if(corr_list[i] > max_corr[1] && dim(samples)[1] > 3) {
    max_corr <- c(corr_list[i], sname, date)
  }
  if(corr_list[i] < min_corr[1] && dim(samples)[1] > 3) {
    min_corr <- c(corr_list[i], sname, date)
  }
  cat("Average correlation between replicates:",corr_list[i],"\n")
}
print(min_corr)
print(max_corr)
plot_data <- data.frame(x=corr_list)
p <- ggplot(plot_data, aes(x)) +
  geom_density() +
  theme_minimal() +
  xlab("correlation between replicates") +
  xlim(c(-1, 1))
p
ggsave("plots/correlation_between_replicates.png", scale=1.5, width=4, height=4, units="in")

measured_cross_correlation <- c()
# what's the correlation between individuals
for(i in 1:(dim(unique_replicates)[1]-1)) {
  for(j in 2:dim(unique_replicates)[1]) {
    if(i != j) {
      # get set of replicates i
      sname.i <- unique_replicates[i]$sname
      date.i <- unique_replicates[i]$collection_date
      idx <- as.character(get_variable(replicate_samples, "sname")) == sname.i
      samples <- prune_samples(idx, replicate_samples)
      idx <- as.character(get_variable(samples, "collection_date")) == date.i
      samples.i <- prune_samples(idx, samples)
      
      # get set of replicates j
      sname.j <- unique_replicates[j]$sname
      date.j <- unique_replicates[j]$collection_date
      idx <- as.character(get_variable(replicate_samples, "sname")) == sname.j
      samples <- prune_samples(idx, replicate_samples)
      idx <- as.character(get_variable(samples, "collection_date")) == date.j
      samples.j <- prune_samples(idx, samples)
      
      # just use the first of the pairs for now
      sample_lr.i <- log_ratios[,colnames(log_ratios) %in% samples.i$sample_id]
      sample_lr.j <- log_ratios[,colnames(log_ratios) %in% samples.j$sample_id]
      total_corr <- 0
      pairs <- 0
      for(k in 1:(dim(sample_lr.i)[2]-1)) {
        for(m in 2:dim(sample_lr.j)[2]) {
          y.t <- as.vector(sample_lr.i[,k])
          y.tt <- sqrt(y.t%*%y.t)
          y.h <- as.vector(sample_lr.j[,m])
          y.hh <- sqrt(y.h%*%y.h)
          total_corr <- total_corr + (y.t%*%y.h)/(y.tt*y.hh)
          pairs <- pairs + 1
        }
      }
      cat("Avg. correl between",sname.i,"-",date.i,"and",sname.j,"-",date.j,":",(total_corr/pairs),"\n")
      measured_cross_correlation[length(measured_cross_correlation)+1] <- total_corr/pairs
    }
  }
}
plot_data <- data.frame(x=measured_cross_correlation)
p <- ggplot(plot_data, aes(x)) +
  geom_density() +
  theme_minimal() +
  xlab("correlation across unrelated replicates") +
  xlim(c(-1, 1))
p
ggsave("plots/correlation_across_replicates.png", scale=1.5, width=4, height=4, units="in")


# plot replicates (proportion total abundance)

replicate_months <- unlist(lapply(md[md$sample_status==2,"collection_date"]@.Data, function(x) substr(x, 6, 7)))
reps_int <- as.numeric(replicate_months)

month_data <- data.frame(month=as.factor(unlist(lapply(reps_int, translate_month))))
levels(month_data$month) <- c("Feb", "Apr", "May", "Jun", "Jul")
p <- ggplot(month_data, aes(x=month)) + geom_bar()
ggsave("plots/replicate_frequency.png", plot=p, scale=1.5, height=4, width=4, units="in")

# see how much batch variables vary across these time points
batch_vars <- c("plate", "extract_dna_conc_ng", "flow_cell_lane", "library_pool")

replicates <- md[md$sample_status==2,]
replicates <- replicates[order(replicates$collection_date),]
replicates <- replicates[order(replicates$sname),]
write.table(t(replicates[,c("sample_id","collection_date","sname","sex","age","season","plate",
                            "extract_dna_conc_ng","library_pool","flow_cell_lane")]),
            file="temp.txt")

# reuse
replicates <- subset_samples(data, sample_status==2)
#glommed_reps_genus <- glom_counts(replicates, level="genus", NArm=FALSE)
glommed_reps_family <- glom_counts(replicates, level="family", NArm=FALSE)

individual <- "NAI"
#indiv_samples <- subset_samples(glommed_reps_genus, sname==individual)
indiv_samples <- subset_samples(glommed_reps_family, sname==individual)

# collapse all lowly expressed guys into one lump
taxa_sums <- colSums(otu_table(indiv_samples)@.Data)
taxa_IDs <- names(taxa_sums)
taxa_sums <- as.vector(taxa_sums)
merge_list <- taxa_IDs[which(taxa_sums < 25)]
merged_samples <- merge_taxa(indiv_samples, merge_list)
# append a rep # to each collection date, just for plotting
fake_dates <- sample_data(merged_samples)$collection_date
for(i in 1:length(fake_dates)) {
  fake_dates[i] <- paste(fake_dates[i],i,sep="_")
}
sample_data(merged_samples)$collection_date <- fake_dates
plot_timecourse(merged_samples, paste(individual,"_replicate_timecourse",sep=""), legend=T, legend_level="genus")
#plot_timecourse(merged_samples, paste(individual,"_replicate_timecourse",sep=""), legend=T, legend_level="family")
