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
  }
}
cat(max_sname,"has max samples:",max_samples,"\n")

max_sname_samples <- subset_samples(samples.filtered, sname==max_sname)
plot_timecourse(max_sname_samples, "DUX_timecourse", legend=TRUE, legend_level="family")

s1.samples <- subset_samples(samples.filtered, season=="Wet")
s1.samples <- subset_samples(s1.samples, sname==max_sname)
s1.num <- nsamples(s1.samples)
s2.samples <- subset_samples(samples.filtered, season=="Dry")
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
adonis(d ~ label, data=labels)

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
biplot(data.pcx, var.axes=T, scale=0, xlab=PC1, ylab=PC2)

PC1.cool.sv <- sequence_vars[12]
PC2.cool.sv <- sequence_vars[28]

cat("SV 12 is:",paste(as.vector(tax_table(samples.filtered)[PC1.cool.sv,]),sep=" "),"\n")
cat("SV 28 is:",paste(as.vector(tax_table(samples.filtered)[PC2.cool.sv,]),sep=" "),"\n")






























