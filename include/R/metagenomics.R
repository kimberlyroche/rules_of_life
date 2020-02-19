base_path <- "/data/mukherjeelab/rulesoflife"

source(file.path(base_path,"include/R/general.R"))

# ====================================================================================================================
# METAGENOMICS FILE HANDLING
# ====================================================================================================================

# read in PiPhillin data from tab-delimited txt file
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are unordered samples
read_metagenomics <- function(metadata, subset=TRUE) {
  # subset=TRUE -- just read in the "Filtered_enzymes.txt" shortlist
  piphillin_dir <- file.path(data_dir,"Piphillin_20190222")
  if(subset) {
    piphillin_file <- "Filtered_enzymes.txt"
    # removed "Enzymes" header at (1,1) in Filtered_enzymes.txt
    whole_dataset <- as.matrix(read.table(file=file.path(piphillin_dir,piphillin_file), sep='\t',
                                           stringsAsFactors=FALSE, header=TRUE, check.names=FALSE))
    samples <- colnames(whole_dataset)
  } else {
    # I chopped this up into 12 files for testing; loop through the 12
    for(i in 1:12) {
      cat("Reading part",i,"\n")
      piphillin_file <- paste("enzymes_pt",i,".txt",sep="")
      data.piphillin <- as.matrix(read.table(file=file.path(piphillin_dir,piphillin_file), sep='\t',
                                             stringsAsFactors=FALSE, header=TRUE, check.names=FALSE))
      enzymes <- data.piphillin[,1]
      data.piphillin <- data.piphillin[,2:dim(data.piphillin)[2]]
      data.piphillin <- apply(data.piphillin, c(1,2), as.numeric)
      rownames(data.piphillin) <- enzymes
      if(is.null(whole_dataset)) {
        whole_dataset <- data.piphillin
      } else {
        whole_dataset <- rbind(whole_dataset, data.piphillin)
      }
    }
    # replace underscores in samples with dashes as in phyloseq object sample IDs
    samples <- colnames(whole_dataset)
    samples <- gsub("_", "-", samples)
    colnames(whole_dataset) <- samples
  }
  # first, exclude samples not present in *both* data sets
  samples.orig <- metadata$sample_id
  samples.both <- intersect(samples, samples.orig)
  
  subset.metadata <- metadata[metadata$sample_id %in% samples.both,]$sample_id
  subset.dataset <- whole_dataset[,samples.both]
  return(subset.dataset)
}

# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are
# samples ordered by collection date (but not sub-ordered by anything)
intersect_order_metagenomics_samples <- function(data.piphillin) {
  sid_collection_date <- metadata[metadata$sample_id %in% samples, c("sample_id", "collection_date")]
  sid_ordered <- sid_collection_date[order(sid_collection_date$collection_date),c("sample_id", "collection_date")]
  ordered.data.piphillin <- data.piphillin[,sid_ordered$sample_id]
  return(ordered.data.piphillin)
}

sname_subset_idx <- function(metadata, sname) {
  idx.subset <- metadata[metadata$sname == sname, c("sample_id","collection_date")]
  idx.subset <- idx.subset[order(idx.subset$collection_date),] # order by collection date
  return(idx.subset)
}

# subset by individual name (sname) and reorder
# this requires 16S metadata to pull sample IDs associated with that host
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are
# samples ordered by collection date and subsetted to individual sname
subset_metagenomics_sname <- function(metagenomics_data, sname, metadata) {
  idx.subset <- sname_subset_idx(metadata, sname)
  sid_ordered <- idx.subset$sample_id
  subset.metagenomics <- metagenomics_data[,colnames(metagenomics_data) %in% idx.subset$sample_id]
  return(subset.metagenomics)
}

# converts and metagnomics matrix to a tidy array of proportions (for input into a timecourse plot)
metagenomics_proportions_tidy <- function(metagenomics_data, sname, metadata) {
  enzymes <- rownames(metagenomics_data)
  idx.subset <- sname_subset_idx(metadata, sname)
  data.prop <- prop.table(metagenomics_data, 2) # column(sample)-wise proportion
  df.prop <- gather_array(data.prop, "proportion", "enzyme", "sample")
  # the enzyme names are super unwieldy, so we'll omit them but we could think about including a
  # a stub as below
  # replace each sample ID with its collection date for readability, sanity checking!
  for(i in 1:dim(df.prop)[1]) {
    # these are uninterpretable if there are a ton of them, so add an average proportion so
    # we can identify the highly expressed guys!
    enz_id <- as.numeric(df.prop$enzyme[i])
    df.prop$enzyme[i] <- enzymes[enz_id]
    sample_placeholder <- as.numeric(df.prop$sample[i])
    sample_label <- colnames(data.prop)[sample_placeholder]
    sample_date <- idx.subset[sample_label]$collection_date
    df.prop$sample[i] <- idx.subset[sample_label]$collection_date
  }
  # fix the order of dates by converting to factors; not sure why this works but ggplot will reorder
  # the samples field otherwise!
  df.prop.ordered <- df.prop[order(df.prop$sample),]
  # df.prop.ordered$sample <- factor(df.prop.ordered$sample, levels=unique(df.prop.ordered$sample))
  # df.prop.ordered$enzyme <- factor(df.prop.ordered$enzyme, levels=unique(df.prop.ordered$enzyme))
  return(df.prop.ordered)
}

# expects a tidy array with columns |enzyme{factor;number}| |sample{factor;date}| |proportion{float}|
plot_timecourse_metagenomics <- function(metagenomics_prop, save_filename="metagenomics_timecourse", gapped=FALSE, legend=FALSE) {
  na.string <- ".N/A"
  # basically no way to make this legend legible, so haven't tested that

  # make sure the dates are in order (redundant, but)
  df2 <- metagenomics_prop[order(metagenomics_prop$sample),]
  if(gapped) {
    # insert empty samples were gaps of 2 weeks or more exist
    gap.days <- 13
    dates_present <- unique(df2$sample)
    for(d in 1:(length(dates_present)-1)) {
      diff <- as.Date(dates_present[d+1]) - as.Date(dates_present[d])
      next.date <- as.Date(dates_present[d])
      attr(diff, "units") <- "days"
      while(diff > gap.days) {
        next.date <- next.date + gap.days
        df2 <- rbind(df2, list(enzyme=na.string, sample=as.character(next.date), proportion=0))
        diff <- as.Date(dates_present[d+1]) - next.date
      }
    }
    img_width <- 15
    df2 <- df2[order(df2$sample),]
  }
  df2$enzyme <- factor(df2$enzyme, levels=unique(df2$enzyme))
  
  categories <- unique(df2$enzyme)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df2$enzyme)))
  if(gapped) {
    img_width <- 15
  } else {
    img_width <- 10
  }
  img_height <- 4
  p <- ggplot(df2, aes(x=sample, y=proportion, fill=enzyme)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8)) +
    theme(axis.text.x = element_text(size=20))
  if(legend) {
    p <- p + theme(legend.position="bottom")
  } else {
    p <- p + theme(legend.position="none")
  }
  ggsave(paste0(save_filename,".png"), plot=p, scale=2, width=img_width, height=img_height, units="in")
}
