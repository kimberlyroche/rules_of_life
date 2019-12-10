source("include/R/general.R")
source("include/R/GP.R") # for load_outcomes()
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
use_tsne <- FALSE
data <- load_and_filter(level=level)

# filter on annotations first
# get group labels
group_labels <- get_group_labels(data)

outcomes <- read.csv(paste0(data_dir,"fitness/IndividualTraits_ForKim.csv"), header=TRUE)
temp <- outcomes[outcomes$sname %in% over_1,]
temp <- temp[!is.na(temp$LRS_livebirths) & temp$LRS_livebirths > 0 &
                   !is.na(temp$LRS_survbirths) & temp$LRS_survbirths > 0 &
                   !is.na(temp$lifetime_rateLiveBirths) & temp$lifetime_rateLiveBirths > 0 &
                   !is.na(temp$lifetime_rateSurvBirths) & temp$lifetime_rateSurvBirths > 0,]
temp <- temp[!is.na(temp$known_lifespan) & temp$known_lifespan > 0,]
temp <- temp[!is.na(temp$age_first_live_birth) & temp$age_first_live_birth > 0,]
cat(paste0("There are ",nrow(temp)," individuals matching fitness annotation selection criteria!\n"))

filtered_snames <- unique(as.vector(temp$sname))

# filter on minimum counts next
min_count_threshold <- 65

subsetted_data <- subset_samples(data, sname %in% filtered_snames)
md <- sample_data(data)
per_indiv_counts <- md %>%
  group_by(sname) %>%
  count()
indiv_list <- per_indiv_counts[per_indiv_counts$n > min_count_threshold,]$sname
samples <- subset_samples(subsetted_data, sname %in% indiv_list)
cat("Number of samples:",phyloseq::nsamples(samples),"\n")
cat("Taxa in each sample:",ntaxa(samples),"\n")

# show replicates by SID, sname, date
md <- sample_data(samples)
temp <- table(md[,"sname"])
temp

# show group breakdown
table(group_labels[indiv_list])

# estimate average retention of samples (after PAM) assuming uniform selection
# 75 starting samples establish a minimum of about 28-30 retained for a given individual
#png("test.png")
#plot(hist(round(2000*(temp/sum(temp)))))
#dev.off()

print(sort(round(2000*(temp/sum(temp))))[1:5])

cat("Number of samples:",phyloseq::nsamples(samples),"\n")

# we'll choose 2000 of these ^^^ samples

#d <- as.data.frame(rep_meta) %>% select(c(sid, sname, collection_date, season))
#d %>%
#  group_by(sname, season) %>%
#  tally()

md <- sample_data(samples)
sample_counts <- otu_table(samples)@.Data
rownames(sample_counts) <- md$sid
colnames(sample_counts) <- NULL

phylo_dist <- TRUE

if(file.exists("PAM_sample_dist.RData")) {
  cat("Loading existing distance matrix...\n")
  load("PAM_sample_dist.RData")
} else {
  cat("Calculating distance matrix...\n")
  if(phylo_dist) {
    sample_dist <- UniFrac(samples, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
  } else {
    # embed (Aitchison distance)
    sample_clr <- driver::clr(sample_counts + pc) # CLR
    sample_dist <- dist(sample_clr)
  }
  save(sample_dist, file=paste0(output_dir,"PAM_sample_dist.RData"))
}
cat("Calculating embedding coordinates...\n")
if(use_tsne) {
  sample_embed <- Rtsne(sample_dist)
} else {
  sample_embed <- cmdscale(sample_dist)
}

# cluster an re-label
k <- 2000
if(file.exists("PAM_sample_cr.RData")) {
  cat("Loading existing PAM clustering...\n")
  load("PAM_sample_cr.RData")
} else {
  cat("Calculating PAM clustering...\n")
  sample_cr <- Cluster_Medoids(as.matrix(sample_dist), k, verbose=TRUE, threads=4)
  save(sample_cr, file=paste0(output_dir,"PAM_sample_cr.RData"))
}
medoids <- logical(nrow(sample_counts))
medoids[sample_cr$medoid_indices] <- TRUE
if(use_tsne) {
  df <- data.frame(x=sample_embed$Y[,1], y=sample_embed$Y[,2], sname=md$sname, Medoid=medoids)
} else {
  df <- data.frame(x=sample_embed[,1], y=sample_embed[,2], sname=md$sname, Medoid=medoids)
}
rownames(df) <- as.character(seq(1, nrow(df)))

save_label <- paste0(plot_dir,"PAM_",level,"_")
if(use_tsne) {
  save_label <- paste0(save_label, "tsne_")
}
save_label <- paste0(save_label, "diagnostic")
p <- ggplot() +
  geom_point(data=df, aes(x=x, y=y, color=sname), size=1)
ggsave(paste0(save_label,"1.png"), plot=p, scale=1, width=12, height=10, units="in", dpi=100)

p <- ggplot(data=df, aes(color=Medoid)) +
  geom_point(data=df[df$Medoid==FALSE,], aes(x=x, y=y), size=1, alpha=0.5) +
  geom_point(data=df[df$Medoid==TRUE,], aes(x=x, y=y), size=1)
ggsave(paste0(save_label,"1.png"), plot=p, scale=1, width=12, height=10, units="in", dpi=100)

per_indiv_medoids <- df %>%
  filter(Medoid == TRUE) %>%
  group_by(sname) %>% count()

# plot a few high-level distributions
p <- ggplot(per_indiv_medoids, aes(n)) +
  geom_histogram()
ggsave(paste0(plot_dir,"PAM_",level,"_diagnostic_indivcount.png"), plot=p, scale=1, width=6, height=6, units="in", dpi=100)

per_indiv_percent <- data.frame(sname=c(), percent=c())
for(sn in per_indiv_medoids$sname) {
  per_indiv_percent <- rbind(per_indiv_percent, data.frame(
  sname=c(sn),
  percent=c(round(per_indiv_medoids[per_indiv_medoids$sname==sn,]$n/per_indiv_counts[per_indiv_counts$sname==sn,]$n, 2))))
}
p <- ggplot(per_indiv_percent, aes(percent)) +
  geom_histogram()
ggsave(paste0(plot_dir,"PAM_",level,"_diagnostic_indivpercent.png"), plot=p, scale=1, width=6, height=6, units="in", dpi=100)

# plot time courses for 5 random selected individuals

random_indiv <- as.character(per_indiv_medoids$sname[sample(nrow(per_indiv_medoids))[1:5]])
for(indiv in random_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(samples, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(plot_dir,"PAM_",level,"_random_",indiv), gapped=FALSE, legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

# plot time courses for MIN and MAX sampled individuals

min_sample_indiv <- as.vector(per_indiv_medoids$sname[which(per_indiv_medoids$n == min(per_indiv_medoids$n))])
cat("Minimally sampled individuals (PAM):\n")
print(min_sample_indiv)

max_sample_indiv <- as.vector(per_indiv_medoids$sname[which(per_indiv_medoids$n == max(per_indiv_medoids$n))])
cat("Maximally sampled individuals (PAM):\n")
print(max_sample_indiv)

for(indiv in min_sample_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(samples, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(plot_dir,"PAM_",level,"_min_",indiv), gapped=FALSE, legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

for(indiv in max_sample_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(samples, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, paste0(plot_dir,"PAM_",level,"_max_",indiv), gapped=FALSE, legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

