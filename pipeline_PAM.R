library(driver)
library(ClusterR)
#library(cluster)
#library(psych)
#library(Rcpp)
library(Rtsne)

source("include.R")

# read in and filter full data set at this phylogenetic level
glom_data <- load_glommed_data(level="family", replicates=TRUE)
data <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2, verbose=TRUE)
replicates <- subset_samples(data, sample_status==2)
cat("Number of replicate samples:",nsamples(replicates),"\n")
cat("Taxa in each replicate:",ntaxa(replicates),"\n")

# show replicates by SID, sname, date
rep_meta <- sample_data(replicates)
table(rep_meta[,"sid"])

d <- as.data.frame(rep_meta) %>% select(c(sid, sname, collection_date, season))
d %>%
  group_by(sname, collection_date, season) %>%
  tally()

rep_counts <- otu_table(replicates)@.Data
rownames(rep_counts) <- rep_meta$sid
colnames(rep_counts) <- NULL

# embed (Aitchison distance)
rep_clr <- driver::clr(rep_counts + 0.5) # CLR
rep_dist <- dist(rep_clr)
res_tsne <- Rtsne(rep_dist)

# cluster an re-label
k <- 30
res_cr <- Cluster_Medoids(as.matrix(rep_dist), k, verbose=TRUE, threads=4)
medoids <- logical(nrow(rep_counts))
medoids[res_cr$medoid_indices] <- TRUE
df <- data.frame(x=res_tsne$Y[,1], y=res_tsne$Y[,2], SID=rownames(res_mds), Medoid=medoids)
p <- ggplot() +
  geom_point(data=df[df$Medoid==FALSE,], aes(x=x, y=y, color=SID), size=1, alpha=0.5) +
  geom_point(data=df[df$Medoid==TRUE,], aes(x=x, y=y, color=SID), size=2)
p

# if I run this on ACA vs. DUI time series data with k=4, are the results as interpretable as
# ACA-wet, ACA-dry, DUI-wet, DUI-dry?

ACA_data <- subset_samples(data, sname=="ACA")
ACA_season <- sample_data(ACA_data)$season
DUI_data <- subset_samples(data, sname=="DUI")
DUI_season <- sample_data(DUI_data)$season
counts <- rbind(otu_table(ACA_data)@.Data, otu_table(DUI_data)@.Data)
sname_labels <- c(rep("ACA", nsamples(ACA_data)), rep("DUI", nsamples(DUI_data)))
season_labels <- c(ACA_season, DUI_season)
labels <- paste(sname_labels, season_labels)
rownames(counts) <- labels
colnames(counts) <- NULL

# embed (Aitchison distance)
rep_clr <- driver::clr(counts + 0.5) # CLR
rep_dist <- dist(rep_clr)
res_tsne <- Rtsne(rep_dist)
df <- data.frame(x=res_tsne$Y[,1], y=res_tsne$Y[,2], label=labels)
p <- ggplot(df) +
  geom_point(aes(x=x, y=y, color=label))
p

# cluster an re-label
k <- 20
res_cr <- Cluster_Medoids(as.matrix(rep_dist), k, verbose=TRUE, threads=4)
medoids <- logical(nrow(counts))
medoids[res_cr$medoid_indices] <- TRUE
df <- data.frame(x=res_tsne$Y[,1], y=res_tsne$Y[,2], label=labels, Medoid=medoids)
p <- ggplot() +
  geom_point(data=df[df$Medoid==FALSE,], aes(x=x, y=y, color=label), size=1, alpha=0.5) +
  geom_point(data=df[df$Medoid==TRUE,], aes(x=x, y=y, color=label), size=2)
p

# visualize medoids on a timecourse
selected_samples <- sample_data(ACA_data)$sid[medoids]
selected_samples <- selected_samples[!is.na(selected_samples)]
plot_timecourse_phyloseq(ACA_data, paste0("ACA_medoids_k",k), gapped=FALSE, legend=FALSE, selected_samples=selected_samples)

selected_samples <- sample_data(DUI_data)$sid[medoids]
selected_samples <- selected_samples[!is.na(selected_samples)]
plot_timecourse_phyloseq(DUI_data, paste0("DUI_medoids_k",k), gapped=FALSE, legend=FALSE, selected_samples=selected_samples)
