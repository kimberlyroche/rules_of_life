glom_data <- load_glommed_data(level="genus", replicates=TRUE)
filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2)
metadata <- read_metadata(filtered)

data.piphillin <- read_metagenomics(metadata)
metadata.metagenomics <- read_metadata_metagenomics(data.piphillin, filtered, metadata)

d1 <- log(c(otu_table(filtered)@.Data) + 0.65)

df <- data.frame(x=d1)
p <- ggplot(df, aes(x)) +
  geom_density() +
  theme_minimal() +
  xlab("log counts")
ggsave(paste("plots/16S_logdensity.png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")

d2 <- log(c(data.piphillin + 0.65))

df <- data.frame(x=d2)
p <- ggplot(df, aes(x)) +
  geom_density() +
  theme_minimal() +
  xlab("log counts")
ggsave(paste("plots/metagenomics_logdensity.png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")


