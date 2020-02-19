# this is the real version of PAM (partitioning around medoids) applied to ABRP samples
# the output of this is a list of K samples spread around the space of variation (defined
# by a distance metric -- either Aitchison or weighted Unifrac) near optimally
#
# output are selected samples and a sanity-check ordination

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))
source(file.path(relative_path,"include/R/GP.R")) # for load_outcomes()
source(file.path(relative_path,"include/R/data_transform.R"))
source(file.path(relative_path,"include/R/visualization.R"))

# additional clustering/ordination stuff
library(ClusterR)
library(Rtsne) # t-SNE if desired

# ======================================================================================
#   filter to individuals with solid fitness annotations
# ======================================================================================

# pass target sample number
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (target sample number)", call.=FALSE)
}
k <- as.numeric(args[1])

# read in data and filter full data set at [level] phylogenetic level
level <- "family"
use_tsne <- FALSE # if TRUE ordinate via t-SNE; if FALSE use PCoA
data <- load_and_filter(level=level)

# get group labels
group_labels <- get_group_labels(data)

outcomes <- read.csv(file.path(relative_path,data_dir,"fitness/IndividualTraits_ForKim.csv"), header=TRUE)
# filter to individuals with > 1 sample
temp <- outcomes[outcomes$sname %in% over_1,]
# filter to individuals with data on live or surviving births
temp <- temp[!is.na(temp$LRS_livebirths) &
                   !is.na(temp$LRS_survbirths) &
                   !is.na(temp$lifetime_rateLiveBirths) &
                   !is.na(temp$lifetime_rateSurvBirths),]
temp <- temp[!is.na(temp$known_lifespan),]
temp <- temp[!is.na(temp$age_first_live_birth),]
cat(paste0("There are ",nrow(temp)," individuals matching fitness annotation selection criteria!\n"))

filtered_snames <- unique(as.vector(temp$sname))
subsetted_data <- subset_samples(data, sname %in% filtered_snames)
n_samples <- phyloseq::nsamples(subsetted_data)

# the minimum sample number is heuristic; we'd like each individual to retain ~35 samples
# how many remain after selection is dynamic
# this seems to be a good choice but we probably want to parameter sweep this
min_count_threshold <- 65

# further filter on DNA concentration >= 3 ng.
subsetted_data <- subset_samples(subsetted_data, extract_dna_conc_ng >= 3)
md <- sample_data(subsetted_data)

per_indiv_counts <- md %>%
  group_by(sname) %>%
  count()
indiv_list <- per_indiv_counts[per_indiv_counts$n > min_count_threshold,]$sname
samples <- subset_samples(subsetted_data, sname %in% indiv_list)
cat("Number of samples:",phyloseq::nsamples(samples),"\n")
cat("Taxa in each sample:",ntaxa(samples),"\n")

# show selected samples grouped by sname
md <- sample_data(samples)
temp <- table(md[,"sname"])
temp

# show selected individuals grouped by grp
table(group_labels[indiv_list])

cat("Number of starting samples:",phyloseq::nsamples(samples),"\n")

# ======================================================================================
#   select samples from these filtered ones
# ======================================================================================

md <- sample_data(samples)
sample_counts <- otu_table(samples)@.Data
rownames(sample_counts) <- md$sid
colnames(sample_counts) <- NULL

phylo_dist <- FALSE # if TRUE use weighted Unifrac; if FALSE use Aitchison distance;
                    # a caution: the root of the tree seems to be randomly initialized
                    # so the sample selection is not fully reproducible

# calculate distances

if(file.exists(file.path(relative_path,output_dir,paste0("PAM_sample_dist_",k,".RData")))) {
  cat("Loading existing distance matrix...\n")
  load(file.path(relative_path,output_dir,paste0("PAM_sample_dist_",k,".RData")))
} else {
  cat("Calculating distance matrix...\n")
  if(phylo_dist) {
    sample_dist <- UniFrac(samples, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
  } else {
    # embed (Aitchison distance)
    sample_clr <- driver::clr(sample_counts + pc) # CLR
    sample_dist <- dist(sample_clr)
  }
  save(sample_dist, file=file.path(relative_path,output_dir,paste0("PAM_sample_dist_",k,".RData")))
}

# calculating embedding for later visualization of included/excluded samples
cat("Calculating embedding coordinates...\n")
if(use_tsne) {
  sample_embed <- Rtsne(sample_dist)
} else {
  sample_embed <- cmdscale(sample_dist)
}

# cluster via PAM
if(file.exists(file.path(relative_path,output_dir,paste0("PAM_clustering_",k,".RData"))))
 {
  cat("Loading existing PAM clustering...\n")
  load(file.path(relative_path,output_dir,paste0("PAM_clustering_",k,".RData")))
} else {
  cat("Calculating PAM clustering...\n")
  sample_cr <- Cluster_Medoids(as.matrix(sample_dist), k, verbose=TRUE, threads=4)
  save(sample_cr, file=file.path(relative_path,output_dir,paste0("PAM_clustering_",k,".RData")))
}

# extract medoids/exemplars and label the embedding with these
medoids <- logical(nrow(sample_counts))
medoids[sample_cr$medoid_indices] <- TRUE
if(use_tsne) {
  df <- data.frame(x=sample_embed$Y[,1], y=sample_embed$Y[,2], sname=md$sname, Medoid=medoids)
} else {
  df <- data.frame(x=sample_embed[,1], y=sample_embed[,2], sname=md$sname, Medoid=medoids)
}
rownames(df) <- as.character(seq(1, nrow(df)))

# ======================================================================================
#   plotting/visualization
# ======================================================================================

save_label <- file.path(relative_path,plot_dir,paste0("PAM_",level,"_"))
if(use_tsne) {
  save_label <- paste0(save_label,"tsne_")
}
save_label <- paste0(save_label,"diagnostic")
p <- ggplot() +
  geom_point(data=df, aes(x=x, y=y, color=sname), size=1)
ggsave(paste0(save_label,"_individual_",k,".png"), plot=p, scale=1, width=12, height=10, units="in", dpi=100)

p <- ggplot(data=df, aes(color=Medoid)) +
  geom_point(data=df[df$Medoid==FALSE,], aes(x=x, y=y), size=2, alpha=0.5) +
  geom_point(data=df[df$Medoid==TRUE,], aes(x=x, y=y), size=2)
ggsave(paste0(save_label,"_medoid_",k,".png"), plot=p, scale=1, width=12, height=10, units="in", dpi=100)

per_indiv_medoids <- df %>%
  filter(Medoid == TRUE) %>%
  group_by(sname) %>% count()

# plot a few high-level distributions
p <- ggplot(per_indiv_medoids, aes(n)) +
  geom_histogram()
ggsave(file.path(relative_path,plot_dir,paste0("PAM_",level,"_diagnostic_indivcount_",k,".png")), plot=p, scale=1, width=6, height=6, units="in", dpi=100)

per_indiv_percent <- data.frame(sname=c(), percent=c())
for(sn in per_indiv_medoids$sname) {
  per_indiv_percent <- rbind(per_indiv_percent, data.frame(
  sname=c(sn),
  percent=c(round(per_indiv_medoids[per_indiv_medoids$sname==sn,]$n/per_indiv_counts[per_indiv_counts$sname==sn,]$n, 2))))
}
p <- ggplot(per_indiv_percent, aes(percent)) +
  geom_histogram()
ggsave(file.path(relative_path,plot_dir,paste0("PAM_",level,"_diagnostic_indivpercent_",k,".png")), plot=p, scale=1, width=6, height=6, units="in", dpi=100)

# plot time courses for 5 random selected individuals
# random_indiv <- as.character(per_indiv_medoids$sname[sample(nrow(per_indiv_medoids))[1:5]])
# for(indiv in random_indiv) {
#   selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
#   selected_samples <- md$sid[selected_indices]
#   indiv_data <- subset_samples(samples, sname==indiv)
#   plot_timecourse_phyloseq(indiv_data, file.path(relative_path,plot_dir,paste0("PAM_",level,"_random_",indiv)), gapped=FALSE,
#                            legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
# }

# plot time courses for individuals with the minimum number of samples selected
# and the maximum number of samples selected

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
  plot_timecourse_phyloseq(indiv_data, file.path(relative_path,plot_dir,paste0("PAM_",level,"_min_",indiv,"_",k)), gapped=FALSE,
                           legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

for(indiv in max_sample_indiv) {
  selected_indices <- as.numeric(rownames(df[df$sname==indiv & df$Medoid==TRUE,]))
  selected_samples <- md$sid[selected_indices]
  indiv_data <- subset_samples(samples, sname==indiv)
  plot_timecourse_phyloseq(indiv_data, file.path(relative_path,plot_dir,paste0("PAM_",level,"_max_",indiv,"_",k)), gapped=FALSE,
                           legend=TRUE, legend_level="family", selected_samples=selected_samples, are_replicates=FALSE)
}

