# quick and dirty - adolescent vs. adult sample correlation
source("include.R")
library(coda.base)

pairwise_corr <- function(samples) {
  corr_vec <- c()
  for(s1 in 1:dim(samples)[1]) {
    for(s2 in 1:dim(samples)[1]) {
      if(s2 > s1) {
        # only do lower triangular
        corr_vec <- c(corr_vec, cor(samples[s1,], samples[s2,]))
      }
    }
  }
  return(corr_vec)
}

glom_data <- load_glommed_data(level="genus", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
# family-agglomerated has replicates; remove these
non_reps <- prune_samples(sample_data(filtered)$sample_status==0, filtered)
metadata <- read_metadata(non_reps)
cat("Filtered down to",nsamples(non_reps),"samples x",ntaxa(non_reps),"taxa\n")

adult_threshold <- 5
min_sample_count <- 5

# which individuals were profiled as adolescents
individuals <- unique(sample_data(non_reps)$sname)
matching_criteria <- c()
for(indiv in individuals) {
  indiv_samples <- prune_samples(sample_data(non_reps)$sname==indiv, non_reps)
  if(min(sample_data(indiv_samples)$age) < adult_threshold && max(sample_data(indiv_samples)$age) > adult_threshold) {
    adolescent_samples <- prune_samples(sample_data(indiv_samples)$age < adult_threshold, indiv_samples)
    per_season <- table(sample_data(adolescent_samples)$season)
    if(dim(per_season) == 2 & per_season['Wet'] > min_sample_count & per_season['Dry'] > min_sample_count) {
      adult_samples <- prune_samples(sample_data(indiv_samples)$age >= adult_threshold, indiv_samples)
      per_season <- table(sample_data(adult_samples)$season)
      if(dim(per_season) == 2 & per_season['Wet'] > min_sample_count & per_season['Dry'] > min_sample_count) {
        cat(indiv,"has sample counts:",nsamples(adolescent_samples),",",nsamples(adult_samples),"\n")
        matching_criteria <- c(matching_criteria, indiv)
      }
    }
  }
}
cat("Total matching criteria:",length(matching_criteria),"\n")

# for comparing correlation

if(FALSE) {
  log_ratios <- apply_ilr(non_reps)
  log_ratios <- scale(log_ratios, center=T, scale=F)
  
  ado <- c()
  adu <- c()
  for(indiv in matching_criteria) {
    sample_info <- metadata[metadata$sname == indiv, c("sample_id", "season", "age")]
    sample_lr <- log_ratios[rownames(log_ratios) %in% sample_info$sample_id,,drop=F]
    adolescent_wet <- log_ratios[rownames(log_ratios) %in% sample_info[(sample_info$season=="Wet" & sample_info$age < adult_threshold),]$sample_id,,drop=F]
    ado <- c(ado, pairwise_corr(adolescent_wet))
    adolescent_dry <- log_ratios[rownames(log_ratios) %in% sample_info[(sample_info$season=="Dry" & sample_info$age < adult_threshold),]$sample_id,,drop=F]
    ado <- c(ado, pairwise_corr(adolescent_dry))
    adult_wet <- log_ratios[rownames(log_ratios) %in% sample_info[(sample_info$season=="Wet" & sample_info$age >= adult_threshold),]$sample_id,,drop=F]
    adu <- c(adu, pairwise_corr(adult_wet))
    adult_dry <- log_ratios[rownames(log_ratios) %in% sample_info[(sample_info$season=="Dry" & sample_info$age >= adult_threshold),]$sample_id,,drop=F]
    adu <- c(adu, pairwise_corr(adult_dry))
  }
  
  df <- data.frame(correlation=ado, age_grp="adolescent")
  df <- rbind(df, data.frame(correlation=adu, age_grp="adult"))
  p <- ggplot(df) +
    geom_density(aes(x=correlation, group=age_grp, fill=age_grp), alpha=0.5) +
    theme_minimal()
  ggsave("plots/adolescent_vs_adult_correlation.png", plot=p, scale=1.5, width=5, height=3, units="in")
}

if(TRUE) {
  pseudocount <- 0.65
  ado <- c()
  adu <- c()
  for(indiv in matching_criteria) {
    cat("Evaluating",indiv,"\n")
    sample_info <- metadata[metadata$sname == indiv, c("sample_id", "season", "age")]
    adolescent_wet <- prune_samples(sample_data(non_reps)$sample_id %in% sample_info[(sample_info$season=="Wet" & sample_info$age < adult_threshold),]$sample_id, non_reps)
    ado <- c(ado, c(dist((otu_table(adolescent_wet)@.Data + pseudocount), method='aitchison')))
    adolescent_dry <- prune_samples(sample_data(non_reps)$sample_id %in% sample_info[(sample_info$season=="Dry" & sample_info$age < adult_threshold),]$sample_id, non_reps)
    ado <- c(ado, c(dist((otu_table(adolescent_dry)@.Data + pseudocount), method='aitchison')))
    adult_wet <- prune_samples(sample_data(non_reps)$sample_id %in% sample_info[(sample_info$season=="Wet" & sample_info$age >= adult_threshold),]$sample_id, non_reps)
    adu <- c(adu, c(dist((otu_table(adult_wet)@.Data + pseudocount), method='aitchison')))
    adult_dry <- prune_samples(sample_data(non_reps)$sample_id %in% sample_info[(sample_info$season=="Dry" & sample_info$age >= adult_threshold),]$sample_id, non_reps)
    adu <- c(adu, c(dist((otu_table(adult_dry)@.Data + pseudocount), method='aitchison')))
  }
  df <- data.frame(distance=ado, age_grp="adolescent")
  df <- rbind(df, data.frame(distance=adu, age_grp="adult"))
  p <- ggplot(df) +
    geom_density(aes(x=distance, group=age_grp, fill=age_grp), alpha=0.5) +
    theme_minimal()
  ggsave("plots/adolescent_vs_adult_Aitchison_dist.png", plot=p, scale=1.5, width=5, height=3, units="in")
}














