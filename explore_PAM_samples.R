source("include/R/general.R")
source("include/R/data_transform.R")
source("include/R/visualization.R")

# additional clustering/ordination stuff
library(ClusterR)
#library(cluster)
#library(psych)
#library(Rcpp)
library(Rtsne)

# read in and filter full data set at this phylogenetic level
level <- "family"
glom_data <- load_glommed_data(level=level, replicates=TRUE)
subsetted_data <- subset_samples(glom_data, sname %in% over_75)
samples <- filter_data(subsetted_data, verbose=TRUE)
cat("Number of replicate samples:",phyloseq::nsamples(samples),"\n")
cat("Taxa in each replicate:",ntaxa(samples),"\n")

# show replicates by SID, sname, date
md <- sample_data(samples)
table(md[,"sname"])

#d <- as.data.frame(rep_meta) %>% select(c(sid, sname, collection_date, season))
#d %>%
#  group_by(sname, season) %>%
#  tally()

sample_counts <- otu_table(samples)@.Data
rownames(sample_counts) <- md$sid
colnames(sample_counts) <- NULL

# embed (Aitchison distance)
sample_clr <- driver::clr(sample_counts + pc) # CLR
if(file.exists("PAM_sample_dist.RData")) {
  cat("Loading existing distance matrix...\n")
  load("PAM_sample_dist.RData")
} else {
  cat("Calculating distance matrix...\n")
  sample_dist <- dist(sample_clr)
  #save(sample_dist, file="PAM_sample_dist.RData")
}
cat("Calculating t-SNE embedding coordinates...\n")
sample_tsne <- Rtsne(sample_dist)

# cluster an re-label
k <- 2000
if(file.exists("PAM_sample_cr.RData")) {
  cat("Loading existing PAM clustering...\n")
  load("PAM_sample_cr.RData")
} else {
  cat("Calculating PAM clustering...\n")
  sample_cr <- Cluster_Medoids(as.matrix(sample_dist), k, verbose=TRUE, threads=4)
  #save(sample_cr, file="PAM_sample_cr.RData")
}
medoids <- logical(nrow(sample_counts))
medoids[sample_cr$medoid_indices] <- TRUE
df <- data.frame(x=sample_tsne$Y[,1], y=sample_tsne$Y[,2], sname=md$sname, Medoid=medoids)

per_indiv_medoids <- df %>%
  filter(Medoid == TRUE) %>%
  group_by(sname) %>% count()

min_sample_indiv <- as.vector(per_indiv_medoids$sname[which(per_indiv_medoids$n == min(per_indiv_medoids$n))])
cat("Minimally sampled individuals (PAM):\n")
print(min_sample_indiv)

max_sample_indiv <- as.vector(per_indiv_medoids$sname[which(per_indiv_medoids$n == max(per_indiv_medoids$n))])
cat("Maximally sampled individuals (PAM):\n")
print(max_sample_indiv)

for(indiv in min_sample_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(subsetted_data, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(plot_dir,"PAM_",level,"_min_",indiv), gapped=FALSE, legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

for(indiv in max_sample_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(subsetted_data, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(plot_dir,"PAM_",level,"_max_",indiv), gapped=FALSE, legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

p <- ggplot() +
  geom_point(data=df[df$Medoid==FALSE,], aes(x=x, y=y, color=sname), size=1, alpha=0.5) +
  geom_point(data=df[df$Medoid==TRUE,], aes(x=x, y=y, color=sname), size=2)
ggsave(paste0(plot_dir,"PAM_",level,"_tsne_diagnostic.png"), plot=p, scale=1, width=10, height=10, units="in", dpi=100)

