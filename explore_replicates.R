source("include/R/general.R")
source("include/R/data_transform.R")
source("include/R/visualization.R")

translate_month <- function(x) {
  if(x == 2) { return("Feb") }
  if(x == 3) { return("Mar") }
  if(x == 4) { return("Apr") }
  if(x == 5) { return("May") }
  if(x == 6) { return("Jun") }
  if(x == 7) { return("Jul") }
}

# ====================================================================================================================
# correlation WITHIN replicate sets
# ====================================================================================================================

level <- "genus"
data <- load_and_filter(level)
replicates <- subset_samples(data, sample_status==2) # replicates only
md <- read_metadata(replicates)
unique_replicates <- unique(md[,c("sname","collection_date")])

# apply log ratio transform
log_ratios <- apply_ilr(replicates)
log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))

# ====================================================================================================================
# correlation WITHIN replicate sets
# ====================================================================================================================

corr_list <- c()
min_corr <- c(Inf, "", "")
max_corr <- c(-Inf, "", "")
for(i in 1:dim(unique_replicates)[1]) {
  sname <- unique_replicates[i]$sname
  date <- unique_replicates[i]$collection_date
  idx <- as.character(get_variable(replicates, "sname")) == sname
  samples <- prune_samples(idx, replicates)
  idx <- as.character(get_variable(samples, "collection_date")) == date
  samples <- prune_samples(idx, samples)
  sample_lr <- log_ratios[,colnames(log_ratios) %in% sample_data(samples)$sample_id]
  total_corr <- 0
  pairs <- 0
  for(j in 1:(dim(sample_lr)[2]-1)) {
    for(k in 2:dim(sample_lr)[2]) {
      if(j != k) {
        #y.t <- as.vector(sample_lr[,j])
        #y.tt <- sqrt(y.t%*%y.t)
        #y.h <- as.vector(sample_lr[,k])
        #y.hh <- sqrt(y.h%*%y.h)
        #cat("Corr:",((y.t%*%y.h)/(y.tt*y.hh)),"\n")
        this_corr <- cor(c(sample_lr[,j]), c(sample_lr[,k]))
        total_corr <- total_corr + this_corr
        pairs <- pairs + 1
      }
    }
  }
  corr_list[i] <- total_corr/pairs
  if(corr_list[i] > max_corr[1] && ntaxa(samples) > 3) {
    max_corr <- c(corr_list[i], sname, date)
  }
  if(corr_list[i] < min_corr[1] && ntaxa(samples)[1] > 3) {
    min_corr <- c(corr_list[i], sname, date)
  }
  cat("Average correlation between replicates for",sname,":",corr_list[i],"\n")
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
ggsave(paste0(plot_dir,"correlation_between_replicates.png"), scale=1.5, width=4, height=4, units="in")

# what is up with lower-correlation replicates?
# check the counts-per-replicates for low-correlation sets of replicates
# e.g. VIB (0.34), DUD (0.52)
lc_reps <- subset_samples(filtered, sname=="VIB")
lc_reps <- subset_samples(lc_reps, sample_status==2)
if(length(unique(sample_data(lc_reps)$collection_date))==1) {
  sids <- sample_data(lc_reps)$"sample_id"
  for(i in 1:length(sids)) {
    total_counts <- sum(otu_table(subset_samples(filtered, sample_id==sids[i]))@.Data)
    cat("Counts rep",i,":",total_counts,"\n")
  }
}

# ====================================================================================================================
# correlation BETWEEN replicate sets
# ====================================================================================================================

measured_cross_correlation <- c()
for(i in 1:(dim(unique_replicates)[1]-1)) {
  for(j in 2:dim(unique_replicates)[1]) {
    if(i != j) {
      # get set of replicates i
      sname.i <- unique_replicates[i]$sname
      date.i <- unique_replicates[i]$collection_date
      idx <- as.character(get_variable(replicates, "sname")) == sname.i
      samples <- prune_samples(idx, replicates)
      idx <- as.character(get_variable(samples, "collection_date")) == date.i
      samples.i <- prune_samples(idx, samples)
      
      # get set of replicates j
      sname.j <- unique_replicates[j]$sname
      date.j <- unique_replicates[j]$collection_date
      idx <- as.character(get_variable(replicates, "sname")) == sname.j
      samples <- prune_samples(idx, replicates)
      idx <- as.character(get_variable(samples, "collection_date")) == date.j
      samples.j <- prune_samples(idx, samples)
      
      # just use the first of the pairs for now
      sample_lr.i <- log_ratios[,colnames(log_ratios) %in% sample_data(samples.i)$sample_id]
      sample_lr.j <- log_ratios[,colnames(log_ratios) %in% sample_data(samples.j)$sample_id]
      total_corr <- 0
      pairs <- 0
      for(k in 1:(dim(sample_lr.i)[2]-1)) {
        for(m in 2:dim(sample_lr.j)[2]) {
          #y.t <- as.vector(sample_lr.i[,k])
          #y.tt <- sqrt(y.t%*%y.t)
          #y.h <- as.vector(sample_lr.j[,m])
          #y.hh <- sqrt(y.h%*%y.h)
          this_corr <- cor(c(sample_lr.i[,k]), c(sample_lr.i[,m]))
          total_corr <- total_corr + this_corr
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
ggsave(paste0(plot_dir,"correlation_across_replicates.png"), scale=1.5, width=4, height=4, units="in")

# ====================================================================================================================
# plot selected replicate sets
# ====================================================================================================================

individual <- "VIB"
indiv_samples <- subset_samples(replicates, sname==individual)

# collapse all lowly expressed guys into one lump
taxa_sums <- colSums(otu_table(indiv_samples)@.Data)
taxa_IDs <- names(taxa_sums)
taxa_sums <- as.vector(taxa_sums)
merge_list <- taxa_IDs[which(taxa_sums < 25)]
merged_samples <- merge_taxa(indiv_samples, merge_list)
mod_dates <- sample_data(merged_samples)$collection_date
for(i in 1:length(fake_dates)) {
  # append a rep # to each collection date, just for plotting
  mod_dates[i] <- paste(fake_dates[i],i,sep="_rep")
}
sample_data(merged_samples)$collection_date <- mod_dates
plot_timecourse_phyloseq(merged_samples, paste0(plot_dir,individual,"_replicate_timecourse"), gapped=F, legend=T, legend_level="family")


