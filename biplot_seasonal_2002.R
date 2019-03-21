library(ggbiplot)
library(coda.base)
library(vegan)

source("include.R")
load("glom_data_genus_reps.RData")

year <- 2002
samples.year <- subset_samples(glom_data, grepl(year, collection_date))
samples.filtered <- filter_data(data=samples.year, verbose=T)
cat(nsamples(samples.filtered),"in",year,"\n")
# grab sequence variant from this output
other_sv <- "AACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGCAGGCGGGAAGACAAGTTGGAAGTGAAAACCATGGGCTCAACCCATGAATTGCTTTCAAAACTGTTTTTCTTGAGTAGTGCAGAGGTAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAGGAACACCAGTGGCGAAGGCGGTCTACTGGGCACCAACTGACGCTGAGGCTCGAAAGCATGGGTAGCAAACAGG"
cat(nsamples(samples.filtered),"\n")
cat(ntaxa(samples.filtered),"\n")

snames <- unique(sample_data(samples.filtered)$sname)
max_samples <- 0
max_sname <- ""
for(s in snames) {
  sname.samples <- subset_samples(samples.filtered, sname==s)
  if(nsamples(sname.samples) > max_samples) {
    max_sname <- s
    max_samples <- nsamples(sname.samples)
    cat(max_sname,"has",max_samples,"samples\n")
  }
}

# remove day-duplicate samples for the purposes of visualization
max_sname_samples <- subset_samples(samples.filtered, sname==max_sname)
metadata <- sample_data(max_sname_samples)
prev_col_date <- "1900-01-01"
keep_ids <- c()
season_labels <- c()
for(s in metadata$sample_id) {
  date <- sample_data(subset_samples(max_sname_samples, sample_id==s))$collection_date
  if(date != prev_col_date) {
    keep_ids[length(keep_ids)+1] <- s
    if(date < "2002-06-01" || date > "2002-11-01") {
      season_labels[length(season_labels)+1] <- 1
    } else {
      season_labels[length(season_labels)+1] <- 0
    }
    prev_col_date <- date
  }
}
max_sname_samples <- subset_samples(samples.filtered, sample_id %in% keep_ids)

plot_timecourse(max_sname_samples, paste(max_sname,"_timecourse",sep=""), legend=TRUE, legend_level="family")

s1.samples <- subset_samples(max_sname_samples, season=="Wet")
s1.samples <- subset_samples(s1.samples, sname==max_sname)
s1.num <- nsamples(s1.samples)
s2.samples <- subset_samples(max_sname_samples, season=="Dry")
s2.samples <- subset_samples(s2.samples, sname==max_sname)
s2.num <- nsamples(s2.samples)

pseudocount <- 0.65

s1.counts <- otu_table(s1.samples)@.Data + pseudocount
s2.counts <- otu_table(s2.samples)@.Data + pseudocount
counts.all <- rbind(s1.counts, s2.counts)
sequence_vars <- colnames(counts.all)
cat("Other is:",which(sequence_vars == other_sv),"\n")

sample_ids <- rownames(counts.all)
colnames(counts.all) <- seq(1,length(sequence_vars))
rownames(counts.all) <- seq(1,length(sample_ids))
d <- dist(counts.all, method='aitchison')

# does season label explain a difference
labels <- data.frame(label=c(rep(0, s1.num), rep(1, s2.num)))
adonis(d ~ label, data=labels, permutations=10000)

# BIPLOT; example from https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot
data.clr <- t(apply(counts.all, 1, function(x){log(x) - mean(log(x))}))
data.pcx <- prcomp(data.clr)
data.mvar <- sum(data.pcx$sdev^2) # sum the total variance
# calculate the PC1 and PC2 variance
PC1 <- paste("PC1: ", round(sum(data.pcx$sdev[1]^2)/data.mvar, 3))
PC2 <- paste("PC2: ", round(sum(data.pcx$sdev[2]^2)/data.mvar, 3))
PC3 <- paste("PC3: ", round(sum(data.pcx$sdev[3]^2)/data.mvar, 3))

PC1
PC2
PC3

# these should separate but if re-run no taxa should repeatedly correlate with separation
season_labels[which(season_labels == 1)] <- "wet"
season_labels[which(season_labels == 0)] <- "dry"
ggbiplot(data.pcx, labels=season_labels, groups=season_labels, ellipse=T)
ggsave("biplot_DUX_2002.png", scale=2, width=5, height=4, units="in")


cool.sv <- sequence_vars[12]
cool.sv <- sequence_vars[33]
cool.sv <- sequence_vars[28]
cool.sv <- sequence_vars[36]
cat("PC1 SV is:",paste(as.vector(tax_table(samples.filtered)[cool.sv,]),sep=" "),"\n")






























