library(cluster)
library(Rcpp)
source("include.R")
sourceCpp("fastCorr.cpp")

filtered <- filter_data()
#sample_data <- subset_samples(filtered, sname %in% c("DUI", "ACA"))
#load("sample_data.RData")
#md <- read_metadata(sample_data)

# sample data is ordered by sname, then collection_date

sample_ilr <- apply_ilr(filtered)
sample_ilr <- t(apply(sample_ilr, 1, function(x) x - mean(x)))

# faster to do this by hand?
corr_mat <- fastCorr(sample_ilr)
png("test_whole_corr.png")
image(corr_mat[1:5000,1:5000])
dev.off()

#diss_mat <- 1-corr_mat

#res <- pam(diss_mat, 10, metric="euclidean", do.swap=TRUE)

# any enrichment for season?
#c1 <- names(which(res$clustering == 1))
#c0 <- names(which(res$clustering != 1))

#contingency <- matrix(0, nrow=2, ncol=2, dimnames=list(c("Cluster1", "NotCluster1"), c("Dry","Wet")))
#contingency[1,1] <- nsamples(md[md$sample_id %in% c1 & md$season=="Dry",]) # cluster 1, dry
#contingency[1,2] <- nsamples(md[md$sample_id %in% c1 & md$season=="Wet",]) # cluster 1, dry
#contingency[2,1] <- nsamples(md[md$sample_id %in% c0 & md$season=="Dry",]) # cluster 1, dry
#contingency[2,2] <- nsamples(md[md$sample_id %in% c0 & md$season=="Wet",]) # cluster 1, dry

#fisher.test(contingency, alternative="greater")
