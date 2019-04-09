library(RColorBrewer)
library(ggplot2)

source("include.R")

# read in PiPhillin data from .tsv file
# requires 16S metadata so we can pare down the samples to those in-common between data sets
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are
# samples ordered by collection date (but not sub-ordered by anything)
read_metagenomics <- function(metadata) {
  piphillin_dir <- "C:/Users/Kim/Documents/rules_of_life/original_data/Piphillin_20190222"
  piphillin_file <- "Filtered_enzymes.txt"
  
  # removed "Enzymes" header at (1,1) in Filtered_enzymes.txt
  data.piphillin <- as.matrix(read.table(file=paste(piphillin_dir,piphillin_file,sep="/"), sep='\t',
                                         stringsAsFactors=FALSE, header=TRUE, check.names=FALSE))
  enzymes <- data.piphillin[,1]
  data.piphillin <- data.piphillin[,2:dim(data.piphillin)[2]]
  data.piphillin <- apply(data.piphillin, c(1,2), as.numeric)
  samples <- colnames(data.piphillin)
  rownames(data.piphillin) <- enzymes
  
  # first, exclude samples not present in *both* data sets
  samples.orig <- metadata$sample_id
  samples.both <- intersect(samples, samples.orig)
  
  subset.metadata <- metadata[metadata$sample_id %in% samples.both,]
  subset.data.piphillin <- data.piphillin[,samples.both]
  
  sid_collection_date <- metadata[metadata$sample_id %in% samples, c("sample_id", "collection_date")]
  sid_ordered <- sid_collection_date[order(sid_collection_date$collection_date),c("sample_id", "collection_date")]
  
  ordered.data.piphillin <- subset.data.piphillin[,sid_ordered$sample_id]
  return(ordered.data.piphillin)
}

# subset by individual name (sname) and reorder
# this requires 16S metadata to pull sample IDs associated with that host
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are
# samples ordered by collection date and subsetted to individual sname
subset_metagenomics_sname <- function(metagenomics_data, sname, metadata) {
  idx.subset <- metadata[metadata$sname == sname, c("sample_id","collection_date")]
  idx.subset <- idx.subset[order(idx.subset$collection_date),] # order by collection date
  sid_ordered <- idx.subset$sample_id
  subset.metagenomics <- metagenomics_data[,colnames(metagenomics_data) %in% idx.subset$sample_id]
  return(subset.metagenomics)
}

# converts and metagnomics matrix to a tidy array of proportions (for input into a timecourse plot)
metagenomics_proportions_tidy <- function(metagenomics_data) {
  data.prop <- prop.table(metagenomics_data, 2)
  df.prop <- gather_array(data.prop, "proportion", "enzyme", "sample")
  # the enzyme names are super unwieldy, so we'll omit them but we could think about including a
  # a stub as below
  # replace each sample ID with its collection date for readability, sanity checking!
  for(i in 1:dim(df.prop)[1]) {
    # df.prop$enzyme[i] <- enzymes[i]
    sample_placeholder <- as.numeric(df.prop$sample[i])
    sample_label <- colnames(data.prop)[sample_placeholder]
    sample_date <- idx.subset[sample_label]$collection_date
    df.prop$sample[i] <- idx.subset[sample_label]$collection_date
  }
  # fix the order of dates by converting to factors; not sure why this works but ggplot will reorder
  # the samples field otherwise!
  df.prop.ordered <- df.prop[order(df.prop$sample),]
  df.prop.ordered$sample <- factor(df.prop.ordered$sample, levels=unique(df.prop.ordered$sample))
  df.prop.ordered$enzyme <- factor(df.prop.ordered$enzyme, levels=unique(df.prop.ordered$enzyme))
  return(df.prop.ordered)
}

# expects a tidy array with columns |enzyme{factor;number}| |sample{factor;date}| |proportion{float}}
plot_metagenomics_timecourse <- function(metagenomics_prop, save_filename="metagenomics_timecourse", legend=FALSE) {
  p <- ggplot(metagenomics_prop, aes(x=sample, y=proportion, fill=enzyme)) + 
    geom_bar(position="fill", stat="identity") +
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8))
  if(legend) {
    p <- p + theme(legend.position="bottom")
  } else {
    p <- p + theme(legend.position="none")
  }
  ggsave(paste(save_filename,".png",sep=""), plot=p, scale=2, width=10, height=4, units="in")
}

data.16S <- read_data()
md.16S <- read_metadata(data.16S)

data.piphillin <- read_metagenomics(md.16S)
subset.piphillin <- subset_metagenomics_sname(data.piphillin, "ACA", md.16S)
prop.piphillin <- metagenomics_proportions_tidy(subset.piphillin)
plot_metagenomics_timecourse(prop.piphillin)

# add gaps?

# labels for enzymes

# next: autocorrelation module for metagenomics

# next: variance component module for metagenomics






