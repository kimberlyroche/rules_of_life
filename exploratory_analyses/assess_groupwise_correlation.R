# the file estimates within-group and between-group correlation as raw numbers to STDOUT

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

# NOTE: refactored and untested

use_individuals <- T # if TRUE we're asking about within-individual correlation vs. between individual correlation
                     # if FALSE we're asking about within- to between-groups
subset <- T # if TRUE just use a subsample (10) of individuals or groups; for estimation
            # if FALSE use all pairs of individuals or all pairs of groups
use_metagenomics <- T # if TRUE we're asking about correlation between Piphillin samples
                      # if FALSE we're asking about correlation between 16S samples

                      # THIS IS UNFINISHED!
                      # THE FOR LOOP BELOW ASSUMES A PHYLOSEQ OBJECT IN ALL ITS SUBSAMPLING

level <- "species"
data <- load_and_filter(level)
metadata <- read_metadata(data)

if(use_metagenomics) {
  # read metagenomics
  data.piphillin <- read_metagenomics(metadata)
  metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, data, metadata)
  data <- data.piphillin
  metadata <- metadata.metagenomics
}

# apply log ratio transform
log_ratios <- apply_ilr(data)
log_ratios <- scale(log_ratios, center=T, scale=F) # samples x taxa or enzymes

if(use_individuals) {
  groups <- unique(md$sname)
} else {
  groups <- unique(md$grp)
}
if(subset) {
  groups <- groups[sample(length(groups))[1:10]]
}
compare_against <- length(groups)-1
within_corr <- c()
between_corr <- c()
for(i in 1:compare_against) {
  # quick and dirty: skip groups for pairwise comparison if each doesn't have
  # some minimum number of samples (e.g. 100)
  if(use_individuals) {
    representation <- nsamples(subset_samples(data, sname==groups[i]))
    min_representation <- 10
  } else {
    representation <- nsamples(subset_samples(data, grp==groups[i]))
    min_representation <- 50
  }
  if(representation > min_representation) {
    for(j in i:length(groups)) {
      if(use_individuals) {
        representation <- nsamples(subset_samples(data, sname==groups[j]))
      } else {
        representation <- nsamples(subset_samples(data, grp==groups[j]))
      }
      if(representation > min_representation) {
        cat("Comparing",groups[i],"to",groups[j],"...\n")
        for(k in 1:min_representation) {
          # sample group A
          if(use_individuals) {
            sample_set.1 <- subset_samples(data, sname==groups[i])
          } else {
            sample_set.1 <- subset_samples(data, grp==groups[i])
          }
          idx <- sample(nsamples(sample_set.1))[1]
          sample_idx.1 <- read_metadata(sample_set.1)$sample_id[idx]
          sample.1 <- subset_samples(sample_set.1, sample_id==sample_idx.1)
          sample_season.1 <- read_metadata(sample.1)$season
          sample_sname.1 <- read_metadata(sample.1)$sname

          # sample group B
          if(use_individuals) {
            sample_set.2 <- subset_samples(data, sname==groups[j])
          } else {
            sample_set.2 <- subset_samples(data, grp==groups[j])
          }
          # match season
          sample_set.2 <- subset_samples(sample_set.2, season == sample_season.1)
          # exclude samples from the same individual
          if(!use_individuals) {
            sample_set.2 <- subset_samples(sample_set.2, sname != sample_sname.1)
          }
          idx <- sample(nsamples(sample_set.2))[1]
          sample_idx.2 <- read_metadata(sample_set.2)$sample_id[idx]
          sample.2 <- subset_samples(sample_set.2, sample_id==sample_idx.2)

          sample_lr.1 <- log_ratios[rownames(log_ratios) %in% sample_idx.1,]
          sample_lr.2 <- log_ratios[rownames(log_ratios) %in% sample_idx.2,]

          # calculate correlation
          corr_est <- cor(sample_lr.1, sample_lr.2)
          if(i == j) {
            within_corr[length(within_corr)+1] <- corr_est
          } else {
            between_corr[length(between_corr)+1] <- corr_est
          }
        }
      }
    }
  }
}

cat("Average within-group correlation:",mean(within_corr),"\n")
cat("Average between-group correlation:",mean(between_corr),"\n")
