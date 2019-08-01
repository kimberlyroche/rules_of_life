library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)
library(scales)
library(coda)
library(Rcpp)
library(RcppEigen)
library(LaplacesDemon) # various sampling distributions
#library(MCMCpack)
library(RColorBrewer)

sourceCpp("fastCorr.cpp")
sourceCpp("dens_optim.cpp")

# these lists were generated manually
# 10 max-sampled individuals
best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")
# individuals with # samples >= 100
over_100 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
              "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO")
# individuals with # samples >= 50
over_50 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
             "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO", "LAD", "LOB", "WAD", "GAB", "LIW",
             "VIN", "TAL", "VEX", "VEI", "ALE", "MBE", "WHE", "WYN", "LOL", "HOL", "NOB", "VOT", "LYE",
             "HON", "DAG", "DUN", "OTI", "LUI", "OFR", "LAZ", "ONY", "VEL", "ELV", "FAX", "ORI", "EAG",
             "ODE", "NIK", "VAP", "WIP", "LOU", "NOO", "EVA", "EXO", "KOR", "NAR", "VOW", "HYM", "PAI",
             "LAS", "VIO", "WEA", "DOU", "LIZ", "WAS", "ZIB", "QUA", "WEN", "WOB", "WOL", "HOK", "LAV",
             "OBI", "POK", "SOR", "KOL", "ISR", "OMO", "SCE", "AFR", "MON", "NIN", "VEB", "ADD", "VOY",
             "DRO", "LOC", "OJU", "OST", "DUB", "LEI", "VAA", "GAN", "HUM", "LUN", "VIV", "BUC", "LAN",
             "LOX", "HAS", "SNA", "WUA", "YAI", "EGO", "ABB", "CRU", "LOF", "WAB", "ZIZ", "COD", "LEX",
             "RAJ", "KIW", "LAO", "LIB", "NJU", "OBR", "OCE", "POW")
# individuals with at least 2 samples
over_1 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
            "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO", "LAD", "LOB", "WAD", "GAB", "LIW",
            "VIN", "TAL", "VEX", "VEI", "ALE", "MBE", "WHE", "WYN", "LOL", "HOL", "NOB", "VOT", "LYE",
            "HON", "DAG", "DUN", "OTI", "LUI", "OFR", "LAZ", "ONY", "VEL", "ELV", "FAX", "ORI", "EAG",
            "ODE", "NIK", "VAP", "WIP", "LOU", "NOO", "EVA", "EXO", "KOR", "NAR", "VOW", "HYM", "PAI",
            "LAS", "VIO", "WEA", "DOU", "LIZ", "WAS", "ZIB", "QUA", "WEN", "WOB", "WOL", "HOK", "LAV",
            "OBI", "POK", "SOR", "KOL", "ISR", "OMO", "SCE", "AFR", "MON", "NIN", "VEB", "ADD", "VOY",
            "DRO", "LOC", "OJU", "OST", "DUB", "LEI", "VAA", "GAN", "HUM", "LUN", "VIV", "BUC", "LAN",
            "LOX", "HAS", "SNA", "WUA", "YAI", "EGO", "ABB", "CRU", "LOF", "WAB", "ZIZ", "COD", "LEX",
            "RAJ", "KIW", "LAO", "LIB", "NJU", "OBR", "OCE", "POW", "IAG", "MLO", "GYP", "LIT", "OPA",
            "COT", "DIP", "LAW", "RHO", "VOR", "AMA", "AYU", "DUR", "FLA", "OAS", "VIB", "CAB", "CHE",
            "HAV", "LUP", "MIC", "YOB", "PIT", "YAN", "LOZ", "TOG", "BEA", "DUD", "GOM", "HIB", "WAG",
            "ETO", "KEL", "NUT", "WES", "IDI", "ISO", "PRU", "YOG", "ZAI", "AZI", "DIC", "EMI", "KAT",
            "LAX", "LIM", "NET", "VAD", "ANE", "CAD", "ECL", "HAM", "OCT", "TAP", "GAS", "NAW", "RWA",
            "CYC", "HES", "LOM", "NOK", "SHY", "CAR", "DJI", "KAG", "LYM", "MOR", "ORN", "RAN", "SCO",
            "CON", "NOS", "WON", "YAR", "BLU", "DYN", "EAS", "IGU", "OOZ", "WIR", "DEA", "DEG", "GAM",
            "KER", "KUT", "NOZ", "SER", "WIV", "APP", "BAT", "HUC", "LAH", "LAT", "LOA", "LUR", "OBL",
            "THU", "ZIN", "DAM", "FAB", "JOB", "KIJ", "OZO", "RIF", "SAG", "WAW", "NOP", "REX", "VES",
            "WIF", "CED", "HEK", "LUX", "NYL", "ROC", "VAN", "WOO", "AMO", "BUT", "DEP", "ELB", "FIG",
            "KRI", "MAH", "YOK", "ZAN", "DAP", "JAG", "LUT", "SAT", "TUS", "VEH", "VIP", "LYC", "NAI",
            "ARS", "HAD", "KRA", "LUD", "NIR", "NUN", "RIK", "UHU", "WAJ", "WEW", "WIZ", "YOI", "BAG",
            "DHO", "FAC", "KEY", "TER", "VUG", "WUM", "GER", "HOF", "IRI", "MYS", "NES", "ROX", "WER",
            "WEU", "ELD", "FIW", "LAA", "LAQ", "SAW", "SUJ", "ATL", "DOV", "EUR", "FAM", "FIN", "VAZ",
            "WYL", "ZAK", "CAY", "DIB", "FUZ", "GYA", "HIV", "LEC", "LUO", "SOK", "VAU", "WAP", "WEI",
            "WOK", "FRU", "GRA", "PEC", "PEM", "TIN", "DED", "HYD", "LYB", "NOJ", "NOM", "WOZ", "YEL",
            "DIG", "GET", "MOZ", "POA", "RHU", "WEZ", "WYC", "APO", "BUL", "ELA", "KIP", "LUS", "SEI",
            "WEG", "WOV", "WRA", "ZOR", "HAP", "KOM", "PLA", "WEY", "WIW", "YUG", "BOL", "DEC", "LOD",
            "LUV", "LUZ", "NAF", "NAG", "NUM", "NYA", "RIG", "SEE", "ALV", "CAV", "ESI", "GEM", "HEL",
            "LEY", "LID", "LOP", "NYU", "POT", "SAC", "WOC", "WOF", "YOL", "ZAM", "ALB", "BIO", "ERN",
            "HIC", "HIZ", "LUF", "NEH", "NON", "VEG", "ASA", "BOC", "DEF", "EBA", "EQU", "ERI", "HEE",
            "LAF", "LES", "RIS", "SEM", "WYA", "YAK", "BIK", "FIL", "GLO", "HOJ", "HOS", "LOS", "MIA",
            "PUY", "RIV", "SAH", "SAS", "ZIM", "ZUR", "AAR", "APH", "CAC", "COP", "DOJ", "DRI", "KAF",
            "KHA", "LEJ", "MOG", "OKO", "ROW", "TAB", "WOI", "ALM", "CHL", "CRI", "EID", "ENV", "HAG",
            "HEA", "ICA", "KOA", "NEA", "NIO", "NUZ", "NYS", "TRU", "WAM", "WED", "WEM", "WEX", "WHO",
            "ABD", "ACI", "BOY", "COU", "DEO", "DRU", "EDG", "ERM", "EXI", "FAK", "GUI", "ILA", "LOJ",
            "LUM", "RAG", "ROD", "SAD", "UTE", "VAQ", "VIA", "WEK", "WOP", "WRY", "YAP", "ACM", "DHA",
            "DIW", "FED", "FER", "GAL", "GUR", "IVY", "LEU", "LUG", "MES", "MEX", "NER", "PIX", "SAU",
            "SQU", "TIV", "UPE", "URU", "YAY", "ZUC", "ACC", "ADR", "ARA", "BOT", "CEN", "DOT", "EDI",
            "EUG", "GLE", "HEX", "HUD", "KNO", "LAL", "LEF", "LEH", "LIX", "NAK", "NTA", "OPU", "ORB",
            "PSY", "RHA", "SPY", "TAY", "TOR", "TRI", "WUN", "WYO", "YAZ", "YIP", "ANC", "BRO", "EKA",
            "GAE", "GEN", "GOV", "GWA", "HAZ", "HIA", "HIN", "HYR", "JIJ", "LAR", "LIV", "LYL", "MEG",
            "NGA", "NOA", "RAZ", "RUB", "RUE", "SUE", "TOF", "USI", "VAT", "WID", "YAB", "YEN", "YUM", "ZAP", "ZAW")

# ====================================================================================================================
# ACCESSORS
# ====================================================================================================================

# read current 16S data into a phyloseq object
read_data <- function(write_sample=FALSE, replicates=TRUE) {
  # rows are samples, columns are taxa
  if(replicates) {
    cat("Using replicates...\n")
    data <- readRDS("original_data/emp_baboon_pool_T_w_techReps.RDS")
  } else {
    cat("NOT using replicates...\n")
    data <- readRDS("original_data/emp_baboon_NewFiltr.RDS")
  }
  # for now, just remove non-Bacterial domain
  data <- subset_taxa(data, domain=="Bacteria")
  if(write_sample) {
    write.table(otu_table(data)[1:10,1:10], file="otu_sample.txt", sep="\t")
    write.table(tax_table(data)[1:10,], file="tax_sample.txt", sep="\t")
  }
  return(data)
}

# pull metadata from phyloseq object
read_metadata <- function(data, write_sample=FALSE) {
  # metadata: rows are samples, "sample-specific variables" are columns
  metadata <- phyloseq::sample_data(data)
  # sort ASC by collection date
  metadata <- metadata[order(metadata$sname,metadata$collection_date),]
  if(write_sample) {
    write.table(sample_data(data)[1:10,1:10], file="samp_sample.txt", sep="\t")
  }
  return(metadata)
}

load_glommed_data <- function(level="species", replicates=TRUE) {
  if(replicates) {
    filename <- paste0("original_data/glom_data_",level,"_reps.RData")
    if(file.exists(filename)) {
      load(filename)
    } else {
      stop(paste0("Agglomerated data file for ",level," (+replicates) does not exist."))
    }
  } else {
    filename <- paste0("original_data/glom_data_",level,".RData")
    if(file.exists(filename)) {
      load(filename)
    } else {
      stop(paste0("Agglomerated data file for ",level," (-replicates) does not exist."))
    }
  }
  return(glom_data)
}

# returns the (sampled) proportion at/above the given threshold
# so estimated percent zeroes is 1 - get_tiny_counts(data, 1)
get_tiny_counts <- function(data, threshold, use_ns=500) {
  if(phyloseq::nsamples(data) < use_ns) {
    use_ns <- phyloseq::nsamples(data)
  }
  sids <- sample(phyloseq::nsamples(data))[1:use_ns]
  return(sum(apply(otu_table(data)@.Data[sids,], c(1,2), function(x) { x >= threshold } ))/(ntaxa(data)*use_ns))
}

# use phyloseq::tax_glom to collapse taxa to level=level
perform_agglomeration <- function(level="genus", replicates=TRUE) {
  data <- read_data(replicates=replicates)
  cat("Zero counts (original):",(1 - get_tiny_counts(data, 1)),"\n")
  if(!replicates) {
    # remove technical replicates
    data <- subset_samples(data, sample_status==0)
  }
  cat("Agglomerating data...\n")
  glom_data <- glom_counts(data, level=level)
  if(replicates) {
    save(glom_data, file=paste("glom_data_",level,"_reps.RData",sep=""))
  } else {
    save(glom_data, file=paste("glom_data_",level,".RData",sep=""))
  }
  cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
}

# filter data below a count threshold in a minimum number of samples into an "Other" category
# e.g. count_threshold=3, sample_threshold=0.2 filters taxa with no more than a 2-count in 80% of samples into an
# "Other" category, labeled by an arbitrary sequence variant in that "Other" category (see the print statement
# below identifying it for reference)
filter_data <- function(data, count_threshold=3, sample_threshold=0.2, collapse_level=NULL, verbose=FALSE) {
  total_counts <- sum(otu_table(data)@.Data)
  # get indices to collapse
  retained_counts <- sum(otu_table(filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= sample_threshold, prune=T))@.Data)
  #collapse_indices <- !as.logical(filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= sample_threshold, prune=F))
  counts <- otu_table(data)@.Data
  collapse_indices <- apply(counts, 2, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) < sample_threshold)
  if(!is.null(collapse_level)) {
    tt <- tax_table(data)@.Data
    collapse_tax_indices <- apply(tt, 1, function(x) is.na(x[collapse_level]))
    collapse_indices <- collapse_indices | collapse_tax_indices
  }
  collapse_taxnames <- names(collapse_indices[which(collapse_indices == TRUE)])
  merged_data <- merge_taxa(data, collapse_taxnames, 1)
  if(verbose) {
    cat("Collapsing",length(collapse_taxnames),"taxa of",ntaxa(data),"\n")
    cat("\tOther category collapsed into:",collapse_taxnames[1],sep=" ","\n")
    collapse_tidx <- which(rownames(tax_table(merged_data)) == collapse_taxnames[1])
    cat("\tIndex of other in tax table:",collapse_tidx,"\n")
    cat("\tTaxonomy of collapsed:",tax_table(merged_data)@.Data[collapse_tidx,],"\n")
    tax_table(merged_data)@.Data[collapse_tidx,] <- rep("Collapsed",7)
    cat("\tCollapsed counts:",(total_counts-retained_counts),"of",total_counts,"(",(total_counts-retained_counts)/total_counts,"total )\n")
    cat("\tPercent zero-count in data set:",(1 - get_tiny_counts(merged_data, 1)),"\n")
  }
  return(merged_data)
}

grp_by_sname <- function(data) {
  snames <- unique(md$sname)
  return(unlist(lapply(snames, function(x) length(unlist(unique(md[md$sname %in% x,"grp"]))))))
}

# there are samples in the 16S ASV data not present in the metagenomics
# remove these from the metadata paired with the metagenomics data set
read_metadata_metagenomics <- function(metagenomics, full_data, full_metadata) {
  sample_ids <- colnames(metagenomics)
  md.sample_ids <- full_metadata$sample_id
  intersection_ids <- intersect(sample_ids, md.sample_ids)
  subset_metadata <- full_metadata[full_metadata$sample_id %in% intersection_ids,]
  return(subset_metadata)
}

# ====================================================================================================================
# DATA TRANSFORMATION, ETC.
# ====================================================================================================================

# recode date as numeric: number of days since earliest sample date
# these are frequently *fractional* because of leap years!!!
date_to_num <- function(date, baseline="2000-04-21") {
  as.numeric(difftime(date, as.Date(baseline), units="days"))
}

filter_data_omit <- function(count_threshold=3, sample_threshold=0.9, data=NULL) {
  if(is.null(data)) {
    load("glom_data_genus.RData")
    #cat("Zero counts (glommed data):",(1 - get_tiny_counts(glom_data, 1)),"\n")
    filtered <- filter_counts(glom_data, count_threshold, sample_threshold)
    #cat("Zero counts (filtered data):",(1 - get_tiny_counts(filtered, 1)),"\n")
  } else {
    filtered <- filter_counts(data, count_threshold, sample_threshold)
  }
  return(filtered)
}

# uses phyloseq::filter_taxa to filter to taxa above a given count_threshold in at least 
# freq_threshold observations
filter_counts <- function(data, count_threshold, freq_threshold) {
  total_counts <- sum(otu_table(data)@.Data)
  filtered_data <- filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= freq_threshold, TRUE)
  total_counts_filtered <- sum(otu_table(filtered_data)@.Data)
  #cat("Total counts:",total_counts,"filtered to:",total_counts_filtered,"\n")
  #cat("Total taxa:",ntaxa(data),"filtered to:",ntaxa(filtered_data),"\n")
  return(filtered_data)
}

# NOTE: shouldn't agglomerate per individual but for testing, we'll try
glom_counts <- function(data, level="genus", NArm=FALSE) {
  cat("Agglomerating data to",level,"level...\n")
  glom_data <- tax_glom(data, taxrank=level, NArm=NArm)
  cat("Agglomerated data has dimension ntaxa=",ntaxa(glom_data),"x nsamples=",phyloseq::nsamples(glom_data),"\n")
  return(glom_data)
}

apply_proportion <- function(data) {
  return(transform_sample_counts(data, function(x) x / sum(x) ))
}

# returns a data.frame
# if passing in a phyloseq object, expects taxa as columns
# if passing in metaganomics data, expects enzymes as rows
apply_ilr <- function(data, pseudocount=0.65) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data # samples (rows) x taxa (columns)
  } else {
    counts <- data # enzymes (rows) x samples (columns)
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  # driver log-ratio transforms expect samples as rows!
  # output: samples (rows) x D-1 taxa (columns)
  return(driver::ilr(counts))
}

# as above; returns D instead of D-1
apply_clr <- function(data, pseudocount=0.65) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data
  } else {
    counts <- data
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  return(driver::clr(counts))
}

# as above
# d is the index of the reference taxon
apply_alr <- function(data, pseudocount=0.65, d=NULL) {
  if("phyloseq" %in% class(data)) {
    counts <- otu_table(data)@.Data
  } else {
    counts <- data
    counts <- t(counts)
  }
  counts <- counts + pseudocount
  return(driver::alr(counts))
}

geom_mean <- function(x) {
  return(prod(x)**(1/length(x)))
}

aitchison_dist <- function(x.i, x.j) {
  sq.sum <- 0
  gm.i <- geom_mean(x.i)
  gm.j <- geom_mean(x.j)
  for(k in 1:length(x.i)) {
    sq.sum <- sq.sum + (log(x.i[k]/gm.i) - log(x.j[k]/gm.j))**2
  }
  return(sqrt(sq.sum))
}

# ====================================================================================================================
# METAGENOMICS FILE HANDLING
# ====================================================================================================================

# read in PiPhillin data from tab-delimited txt file
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are unordered samples
read_metagenomics <- function(metadata, subset=TRUE) {
  # subset=TRUE -- just read in the "Filtered_enzymes.txt" shortlist
  piphillin_dir <- "original_data/Piphillin_20190222"
  if(subset) {
    piphillin_file <- "Filtered_enzymes.txt"
    # removed "Enzymes" header at (1,1) in Filtered_enzymes.txt
    whole_dataset <- as.matrix(read.table(file=paste(piphillin_dir,piphillin_file,sep="/"), sep='\t',
                                           stringsAsFactors=FALSE, header=TRUE, check.names=FALSE))
    samples <- colnames(whole_dataset)
  } else {
    # I chopped this up into 12 files for testing; loop through the 12
    for(i in 1:12) {
      cat("Reading part",i,"\n")
      piphillin_file <- paste("enzymes_pt",i,".txt",sep="")
      data.piphillin <- as.matrix(read.table(file=paste(piphillin_dir,piphillin_file,sep="/"), sep='\t',
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
  ggsave(paste("plots/", save_filename,".png",sep=""), plot=p, scale=2, width=img_width, height=img_height, units="in")
}

# ====================================================================================================================
# VISUALIZATION -- HISTOGRAMS
# ====================================================================================================================

histogram_abundances <- function(data, filename="histogram") {
  df <- gather_array(data)
  p <- ggplot(df) +
    geom_histogram(aes(x=var), binwidth=0.25) +
    theme_minimal() +
    xlab("log ratio abundance")
  ggsave(paste(filename,".png",sep=""), plot=p)
}

histogram_indiv_samples <- function(data) {
  nt <- ntaxa(data)
  per_indiv <- psmelt(data) %>% group_by(sname) %>% summarise(n = n()/nt)
  p <- ggplot(per_indiv) +
    geom_histogram(aes(x=n), binwidth=5) +
    theme_minimal() +
    xlab("per-individual samples")
  ggsave("histogram_per_individual_samples.png", plot=p)
}

histogram_sample_density <- function(data, units="weeks") {
  per_indiv <- psmelt(data)
  snames <- unique(per_indiv$sname)
  differences <- c()
  min_diff <- Inf
  max_diff <- -Inf
  for(s in snames) {
    indiv_sname <- s
    df <- per_indiv[per_indiv$sname==indiv_sname,c("sample_id", "collection_date")]
    df <- unique(df)
    df <- df[order(df$collection_date),"collection_date"]
    if(length(df) > 1) {
      cat("Getting sample differences for individual",s,"...\n")
      indiv_differences <- as.numeric(length(df)-1)
      for(i in 1:(length(df)-1)) {
        d1 <- strptime(df[i], format="%Y-%m-%d")
        d2 <- strptime(df[i+1], format="%Y-%m-%d")
        diff <- d2-d1
        units(diff) <- units
        indiv_differences[i] <- as.numeric(diff)
        if(indiv_differences[i] > max_diff) {
          max_diff <- indiv_differences[i]
          cat("New MAX diff",max_diff,"from individual",s,"\n")
          print(d1)
          print(d2)
        }
        if(indiv_differences[i] < min_diff) {
          min_diff <- indiv_differences[i]
          cat("New MIN diff",min_diff,"from individual",s,"\n")
          print(d1)
          print(d2)
        }
      }
      differences <- c(differences, indiv_differences)
    }
  }
  diff_df <- as.data.frame(x=differences)
  p <- ggplot(diff_df) +
          geom_histogram(aes(x=differences), binwidth=4) +
          theme_minimal() +
          xlab(paste("sample distance (",units,")",sep=""))
  if(units == "weeks") {
    p <- p + xlim(0, 100)
  }
  ggsave("histogram_sample_distance.png", plot=p)
}

# takes a list of correlations or covariances
histogram_corr <- function(data, filename, cov=FALSE) {
  df <- as.data.frame(data)
  colnames(df) <- c("corr")
  p <- ggplot(df) +
         geom_density(aes(x=corr)) +
         theme_minimal()
  if(cov) {
    p <- p + xlab("covariance")
  } else {
    p <- p + xlim(-1, 1) + xlab("correlation")
  }
  ggsave(paste(filename,".png",sep=""), plot=p)
}

# ====================================================================================================================
# VISUALIZATION -- COVARIANCE MATRICES/HEATMAPS
# ====================================================================================================================

plot_corr_matrix <- function(data, filename, cov=FALSE) {
  cor_mat <- fastCorr(t(data))
  #cov_mat <- cov(t(data))
  #cor_mat <- cov2cor(cov_mat)
  df <- gather_array(cor_mat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    scale_fill_gradientn(limits = c(-1,1),
    colours=c("navyblue", "darkmagenta", "darkorange1")) +
    xlab("samples") +
    ylab("samples") +
    theme_minimal()
  if(cov) {
    p <- p + guides(fill=guide_legend(title="covariance"))
  } else {
    p <- p + guides(fill=guide_legend(title="correlation"))
  }
  #ggsave(paste(filename,".pdf",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
  ggsave(paste(filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

visualize_groupwise_covariance <- function(data, md, group, sample=1000000) {

  if(group == "grp") {
    groups <- unique(md$grp); partition_obj <- md$grp; use_no <- sample
  } else if(group == "season") {
    groups <- unique(md$season); partition_obj <- md$season; use_no <- sample
  } else if(group == "matgrp") {
    groups <- unique(md$matgrp); partition_obj <- md$matgrp; use_no <- sample
  } else if(group == "sname") {
    # subsample the individuals
    #groups <- best_sampled; partition_obj <- md$sname; use_no <- sample
    groups <- unique(md$sname); partition_obj <- md$sname; use_no <- sample
  } else if(group == "sex") {
    groups <- unique(md$sex); partition_obj <- md$sex; use_no <- sample
  } else if(group == "age") {
    groups <- c(4.5, 19.0, 30.0); partition_obj <- md$age; use_no <- sample
  } else if(group == "plate") {
    # subsample the plates
    groups <- sample(unique(md$plate))[1:20]; partition_obj <- md$plate; use_no <- sample
  } else if(group == "flow_cell_lane") {
    groups <- unique(md$flow_cell_lane); partition_obj <- md$flow_cell_lane; use_no <- sample
  } else if(group == "library_pool") {
    groups <- unique(md$library_pool); partition_obj <- md$library_pool; use_no <- sample
  } else if(group == "extract_dna_conc_ng") {
    groups <- c(2.5, 5, 7.5, 10, 15, 20, 25, 100); partition_obj <- md$extract_dna_conc_ng; use_no <- sample
  }
  
  # list corresponding to elements of groups with sample IDs in each group's bucket
  partition_idx <- list()
  if(group == "age" || group == "extract_dna_conc_ng") {
    for(g in 1:length(groups)) {
      gmin <- 0
      gmax <- groups[g]
      if(g > 1) {
        gmin <- groups[g-1]
      }
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj > gmin & partition_obj <= gmax,
                                                    "sample_id"][[1]]
    }
  } else {
    for(g in groups) {
      partition_idx[[(length(partition_idx)+1)]] <- md[partition_obj==g, "sample_id"][[1]]
    }
  }

  lr <- apply_ilr(data) # after this data is samples x taxa or enzymes
  lr <- scale(lr, center=T, scale=F)
  partition_samples <- list()
  for(i in 1:length(groups)) {
    if("phyloseq" %in% class(data)) {
      partition_samples[[i]] <- lr[partition_idx[[i]],,drop=FALSE]
    } else {
      partition_samples[[i]] <- lr[rownames(lr) %in% partition_idx[[i]],,drop=FALSE]
    }
  }
  # test clean partition via unique(md[md$sample_id %in% partition_idx[[1]], "grp"]) etc.
  
  stacked_lr <- NULL
  samples <- list()
  for(i in 1:length(groups)) {
    total_samples <- dim(partition_samples[[i]])[1]
    use_ids <- sample(total_samples)
    use_upper <- use_no
    #if(total_samples < use_upper) {
    #  use_upper <- total_samples
    #}
    # omit groups with fewer than the requested number of samples
    if(total_samples >= use_upper) {
      cat(groups[i],"has",total_samples,"samples; using",use_upper,"\n")
      samples[[i]] <- partition_samples[[i]][use_ids[1:use_upper],]
      if(is.null(stacked_lr)) {
        stacked_lr <- samples[[i]]
      } else {
        #stacked_lr <- cbind(stacked_lr, samples[[i]])
        stacked_lr <- rbind(stacked_lr, samples[[i]])
      }
    }
  }
  cat("Sample matrix is",dim(stacked_lr)[1],"by",dim(stacked_lr)[2],"\n")
  if("phyloseq" %in% class(data)) {
    save_filename <- paste("plots/", group, "_cov_matrix", sep="")
  } else {
    save_filename <- paste("plots/", group, "_cov_matrix_metagenomics", sep="")
  }
  save_filename <- 
  plot_corr_matrix(stacked_lr, save_filename)
}

# ====================================================================================================================
# VISUALIZATION -- OTHER
# ====================================================================================================================

# for all taxa, plots the percent at/above the given threshold, ordered descending
plot_percent_threshold <- function(data, threshold=3, save_filename) {
  cat("Plotting percent counts >=",threshold,"in all ASVs...\n")
  count_table <- otu_table(data)
  min_counts <- apply(count_table, 2, function(x) { sum(x >= threshold)/phyloseq::nsamples(data) })
  min_counts <- stack(sort(min_counts, decreasing=TRUE))
  p <- min_counts %>% ggplot(aes(x=seq(1,ntaxa(data)), y=values)) +
    geom_point(size=1) +
    theme_minimal() +
    xlab("ASV no.") +
    ylab(paste("Percent counts >= ",threshold,sep=""))
  ggsave(save_filename, plot=p, scale=1.5, width=5, height=3, units="in")
}

# plot proportional change over time
# note: palette size would probably have to increase to accomodate legend_level finer than family!
plot_timecourse_phyloseq <- function(data, save_filename, gapped=FALSE,
                                     legend=TRUE, legend_level="family", selected_samples=NULL) {
  n <- phyloseq::nsamples(data)
  p <- apply_proportion(data)
  df <- psmelt(p)
  df2 <- bind_cols(list(OTU=df$OTU, Sample=df$Sample, Abundance=df$Abundance, SID=df$sid))
  if(!is.null(selected_samples)) {
    df2$alpha <- 0.25
    df2[df2$SID %in% selected_samples,]$alpha <- 1
  }

  # for gapped plots, this is the OTU ID placeholder we'll use for dummy data
  na.string <- ".N/A"
  
  # replace Sample ID's with their dates for readability
  metadata <- sample_data(data)
  for(i in 1:dim(df2)[1]) {
    df2$Sample[i] <- paste(metadata[metadata$sample_id==df2$Sample[i],"collection_date"][[1]],
                           metadata[metadata$sample_id==df2$Sample[i],"sid"][[1]])
  }

  # make sure the dates are in order and fix the order by converting to factors
  df3 <- df2[order(df2$Sample),]
  if(gapped) {
    # insert empty samples were gaps of 2 weeks or more exist
    gap.days <- 13
    dates_present <- unique(df3$Sample)
    for(d in 1:(length(dates_present)-1)) {
      diff <- as.Date(dates_present[d+1]) - as.Date(dates_present[d])
      next.date <- as.Date(dates_present[d])
      attr(diff, "units") <- "days"
      while(diff > gap.days) {
        next.date <- next.date + gap.days
        df3 <- rbind(df3, list(OTU=na.string, Sample=as.character(next.date), Abundance=1))
        diff <- as.Date(dates_present[d+1]) - next.date
      }
    }
  }
  df3 <- df3[order(df3$Sample),]

  df3$Sample <- factor(df3$Sample, levels=unique(df3$Sample))
  level_idx <- 7
  if(legend_level == "kingdom") { level_idx = 1; }
  if(legend_level == "phylum") { level_idx = 2; }
  if(legend_level == "class") { level_idx = 3; }
  if(legend_level == "order") { level_idx = 4; }
  if(legend_level == "family") { level_idx = 5; }
  if(legend_level == "genus") { level_idx = 6; }
  # replace ASV sequences with their (abbrev.) taxonomy for readability
  for(i in 1:dim(df3)[1]) {
    # show labels as order/family/genus
    # species is NA for all
    if(df3$OTU[i] != na.string) {
      df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],level_idx]),collapse="/")
    }
  }

  categories <- unique(df3$OTU)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df3$OTU)))
  if(gapped) {
    coul[1] <- "#DDDDDD"
    img_width <- 15
  } else {
    img_width <- 10
  }
  if(!is.null(selected_samples)) {
    p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU, alpha=alpha)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8)) +
      theme(legend.position="none")
  } else {
    p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8))
    if(legend) {
      p <- p + theme(legend.position="bottom")
    } else {
      p <- p + theme(legend.position="none")
    }
  }
  if(n < 20) {
    # these are likely replicates
    img_width <- 5
  }
  ggsave(paste("plots/",save_filename,".png",sep=""), plot=p, scale=1.5, dpi=100, width=img_width, height=3, units="in")
}

# this is just a wrapper to call plot_timecourse_phyloseq on a representative sample of individuals
# baboons is a list of snames, e.g. c("ACA", "DUI", "CAI", "COB", "DAS")
perform_mult_timecourse <- function(data, baboons, gapped=FALSE, legend=FALSE) {
  for(b in baboons) {
    b <<- b
    cat("Plotting timecourse for",b,"...\n")
    B_counts <- subset_samples(data, sname==b)
    B_log_ratios <- apply_ilr(B_counts)
    if(gapped) {
      save_filename <- paste("plots/",b[1],"_gapped_timecourse",sep="")
    } else {
      save_filename <- paste("plots/",b[1],"_timecourse",sep="")
    }
    plot_timecourse_phyloseq(B_counts, save_filename=save_filename, gapped=gapped, legend=legend)
  }
}

calc_autocorrelation <- function(data,
                                 metadata,
                                 resample=FALSE,
                                 lag.max=26,
                                 date_diff_units="weeks",
                                 use_lr="clr",
                                 alr_ref=NULL,
                                 resample_rate=0.2) {
  if(date_diff_units != "weeks" && date_diff_units != "months" && date_diff_units != "seasons") {
    date_diff_units <- "weeks"
  }

  individuals <- unique(metadata$sname)
  #individuals <- individuals[1:20] # for testing

  # only if using 
  season_boundaries <- c(200005, 200010, 200105, 200110, 200205, 200210, 200305, 200310,
                       200405, 200410, 200505, 200510, 200605, 200610, 200705, 200710,
                       200805, 200810, 200905, 200910, 201005, 201010, 201105, 201110,
                       201205, 201210, 201305, 201310)

  rounds <- 1
  if(resample) {
    rounds <- 100
    #rounds <- 10 # for testing
  }
  lags <- matrix(0, nrow=lag.max, ncol=rounds)

  # TBD: try various CoDa proportionality measures
  if(use_lr == "clr") {
    log_ratios <- apply_clr(data)
  } else if(use_lr == "alr") {
    log_ratios <- apply_alr(data, d=alr_ref)
  } else {
    log_ratios <- apply_ilr(data)
  }
  log_ratios <- scale(log_ratios, center=TRUE, scale=FALSE)

  for(r in 1:rounds) {
    if(resample) {
      cat("Resampling iteration",r,"\n")
    }
    lag.sums <- numeric(lag.max)
    lag.measured <- numeric(lag.max)
    for(indiv in individuals) {
      cat("Individual",indiv,"...\n")
      # pick an individual
      # this weird syntactic hack seems to be necessary for subset_samples?
      # apparently the thing you're filtering against must be globally available
      indiv <<- indiv
      # pull sample IDs associated with this individual
      sample_info <- metadata[metadata$sname == indiv, c("sample_id", "collection_date")]
      # appear already ordered but let's be paranoid
      sample_info <- sample_info[order(sample_info$collection_date),]
      # we need individuals with at least 2 samples
      if(length(intersect(rownames(log_ratios), sample_info$sample_id)) > 1) {
        sample_lr <- log_ratios[rownames(log_ratios) %in% sample_info$sample_id,,drop=F]
        indiv_sample_no <- dim(sample_lr)[1]
        # randomly censor some of the samples
        do_sample <- rep(1, indiv_sample_no)
        if(resample) {
          do_sample <- rbinom(indiv_sample_no, 1, resample_rate)
        }
        # replace these sample IDs with collection dates, that's what we really want
        # order should be maintained here
        rownames(sample_lr) <- sample_info[sample_info$sample_id %in% rownames(sample_lr), "collection_date"]$collection_date
  
        # get distances between adjacent timepoints in units of {date_diff_units}
        d1 <- as.Date(sample_info$collection_date[1])
        time_diff <- matrix(0, nrow=indiv_sample_no, ncol=indiv_sample_no)
        if(indiv_sample_no > 1) {
          for(d1_idx in 1:(indiv_sample_no-1)) {
            for(d2_idx in (d1_idx+1):indiv_sample_no) {
              d1 <- as.Date(sample_info$collection_date[d1_idx])
              d2 <- as.Date(sample_info$collection_date[d2_idx])
              time_diff[d1_idx,d2_idx] <- as.numeric(difftime(d2, d1, units="weeks"))
              if(date_diff_units == "weeks") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx])
              } else if(date_diff_units == "months") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx]/4)
              } else if(date_diff_units == "seasons") {
                d1_ym <- as.numeric(format(d1, "%Y%m"))
                d2_ym <- as.numeric(format(d2, "%Y%m"))
                # distance is number of season boundaries between the samples, indicated in the vector above
                time_diff[d1_idx,d2_idx] <- sum(season_boundaries < d2_ym) - sum(season_boundaries < d1_ym) + 1
              }
            }
          }
        }
        for(lag in 1:lag.max) {
          # find pairs with this lag
          idx <- which(time_diff==lag, arr.ind=T)
          if(dim(idx)[1] >= 1) {
            for(t in 1:dim(idx)[1]) {
              d1_idx <- idx[t,1]
              d2_idx <- idx[t,2]
              if(do_sample[d1_idx] && do_sample[d2_idx]) {
                # given centering, cosine angle == Pearson's correlation
                lag.sums[lag] <- lag.sums[lag] + cor(sample_lr[d1_idx,], sample_lr[d2_idx,])
                lag.measured[lag] <- lag.measured[lag] + 1
              }
            }
          }
          if(lag.sums[lag] > 0 & lag.measured[lag] > 0) {
            lags[lag,r] <- lag.sums[lag]/lag.measured[lag]
          }
        }
      }
    }
    if(rounds == 1) {
      # print out autocorrelations
      for(lag in 1:lag.max) {
        cat("Lag",lag,"=",lags[lag,r],"\n")
      }
      cat("\n")
    }
  }

  if(resample) {
    # get 90% confidence intervals, summarize
    lower_CI <- as.numeric(lag.max)
    upper_CI <- as.numeric(lag.max)
    avg_ac <- as.numeric(lag.max)
    for(lag in 1:lag.max) {
      ac_measured <- sort(lags[lag,])
      lower_CI[lag] <- ac_measured[round(rounds*0.1)]
      upper_CI[lag] <- ac_measured[round(rounds*0.9)]
      avg_ac[lag] <- mean(ac_measured)
    }
    return(data.frame(cbind(lag=seq(1,lag.max), ac=avg_ac, lower=lower_CI, upper=upper_CI)))
  } else {
    return(as.data.frame(cbind(lag=seq(1,lag.max), ac=as.vector(lags))))
  }
}

plot_mean_autocorrelation <- function(lags, filename="autocorrelation", width=6, height=3) {
  lag.max <- max(lags$lag)
  p <- ggplot(lags, aes(x=lag, y=ac)) +
    geom_line(linetype = "dashed") +
    scale_x_continuous(breaks=seq(0,lag.max,1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  ggsave(paste(filename,".png",sep=""), plot=p, dpi=100, scale=1.5, width=width, height=height, units="in")
}

plot_bounded_autocorrelation <- function(lags, filename="autocorrelation", width=6, height=3) {
  p <- ggplot(lags, aes(x=lag, y=ac)) +
    geom_line(linetype="solid") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    scale_x_continuous(breaks=seq(0,dim(lags)[1],1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  ggsave(paste(filename,".png",sep=""), plot=p, dpi=100, scale=1.5, width=width, height=height, units="in")
}


# ====================================================================================================================
# ESTIMATION OF VARIANCE COMPONENTS
# ====================================================================================================================

# plots crude variance components (per-individual) for diagnostic purposes
plot_cov <- function(datamat, filename) {
  df <- gather_array(datamat)
  p <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("samples") +
    ylab("samples") +
    theme_minimal() +
    guides(fill=guide_legend(title="covariance"))
  ggsave(paste("plots/",filename,".png",sep=""), plot=p, scale=1.5, width=4, height=3, units="in")
}

# pass in zero-filtered data
estimate_variance_components <- function(data, metadata, optim_it=1, use_individuals=Inf, include_residual=False) {
  ilr_data <- NULL
  week_kernel <- NULL
  season_vector <- NULL
  age_vector <- NULL
  group_vector <- NULL
  indiv_vector <- NULL
  plate_vector <- NULL
  conc_vector <- NULL
  
  within_group_corr <- 0.1
  between_group_corr <- 0

  # build kernels for factors of interest

  baboons <- over_50
  if(use_individuals < Inf) {
    baboons <- baboons[sample(length(baboons))[1:use_individuals]] # subsample for debugging
  }
  for(b in 1:length(baboons)) {
    cat("Building kernel for",baboons[b],"\n")
    # get samples for this individual only
    if("phyloseq" %in% class(data)) {
      # phyloseq::subset_samples is causing all kinds of issues, so let's subset on sname manually!
      # basically you can't use subset_samples inside functions or loops
      # see: https://github.com/joey711/phyloseq/issues/752
      remove_idx <- as.character(get_variable(filtered, "sname")) == baboons[b]
      baboon_counts <- prune_samples(remove_idx, filtered)
      # subset metadata
      baboon_metadata <- read_metadata(baboon_counts)
    } else {
      baboon_counts <- subset_metagenomics_sname(data, baboons[b], metadata)
      # column names are sample IDs
      sample_ids <- colnames(baboon_counts)
      baboon_metadata <- metadata[metadata$sample_id %in% sample_ids,]
    }
    baboon_log_ratios <- t(apply_ilr(baboon_counts))
    # we'll ultimately pass ilr_data to optim ax taxa or enzymes (rows) x samples (columns)
    # so build it up in that orientation
    if(is.null(ilr_data)) {
      ilr_data <- baboon_log_ratios
    } else {
      # ilr_data is samples as rows, taxa or enzymes as columns
      ilr_data <- cbind(ilr_data, baboon_log_ratios)
    }

    if("phyloseq" %in% class(data)) {
      nt <- ntaxa(baboon_counts)
      ns <- phyloseq::nsamples(baboon_counts)
    } else {
      nt <- dim(baboon_log_ratios)[1] + 1 # no. taxa = no. enzymes; this is pre-LR transform
      ns <- dim(baboon_log_ratios)[2]
    }

    rho <- 1/2
    subset_week_kernel <- matrix(0, nrow=ns, ncol=ns)
    for(i in 1:ns) {
      for(j in 1:ns) {
        if(i == j) {
          subset_week_kernel[i,j] <- 1
        } else {
          d1 <- as.Date(baboon_metadata$collection_date[i])
          d2 <- as.Date(baboon_metadata$collection_date[j])
          # get number of 'weeks difference'
          week_d <- abs(round((d2-d1)/7))[[1]] + 1
          # distance of 1 is < 7 days
          # distance of 2 is < 14 days, etc.
          subset_week_kernel[i,j] <- rho^week_d + 0.1
        }
      }
    }

    if(is.null(week_kernel)) {
      week_kernel <- subset_week_kernel
    } else {
      prev_r <- dim(week_kernel)[1]
      prev_c <- dim(week_kernel)[2]
      week_kernel <- cbind(week_kernel, matrix(0, nrow=dim(week_kernel)[1], ncol=ns))
      week_kernel <- rbind(week_kernel, matrix(0, nrow=ns, ncol=dim(week_kernel)[2]))
      new_r <- dim(week_kernel)[1]
      new_c <- dim(week_kernel)[2]
      week_kernel[(prev_r+1):new_r,(prev_c+1):new_c] <- subset_week_kernel
    }

    subset_season_vector <- character(ns)
    subset_age_vector <- numeric(ns)
    subset_group_vector <- numeric(ns)
    subset_indiv_vector <- numeric(ns)
    subset_plate_vector <- numeric(ns)
    subset_conc_vector <- numeric(ns)
    for(i in 1:ns) {
      if(baboon_metadata$season[i] == "Wet") {
        subset_season_vector[i] <- "W"
      } else {
        subset_season_vector[i] <- "D"
      }
      if(baboon_metadata$age[i] < 4.5) {
        subset_age_vector[i] <- "J"
      } else if(baboon_metadata$age[i] > 19) {
        subset_age_vector[i] <- "G"
      } else {
        subset_age_vector[i] <- "A"
      }
      subset_group_vector[i] <- baboon_metadata$grp[i]
      subset_indiv_vector[i] <- baboon_metadata$sname[i]
      subset_plate_vector[i] <- baboon_metadata$plate[i]
      if(baboon_metadata$extract_dna_conc_ng[i] <= 2.5) {
        subset_conc_vector[i] <- 1
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 5.0) {
        subset_conc_vector[i] <- 2
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 7.5) {
        subset_conc_vector[i] <- 3
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 10.0) {
        subset_conc_vector[i] <- 4
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 15.0) {
        subset_conc_vector[i] <- 5
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 20.0) {
        subset_conc_vector[i] <- 6
      } else if(baboon_metadata$extract_dna_conc_ng[i] <= 25.0) {
        subset_conc_vector[i] <- 7
      } else {
        subset_conc_vector[i] <- 8
      }
    }
    if(is.null(season_vector)) {
      season_vector <- subset_season_vector
    } else {
      season_vector <- c(season_vector, subset_season_vector)
    }
    if(is.null(age_vector)) {
      age_vector <- subset_age_vector
    } else {
      age_vector <- c(age_vector, subset_age_vector)
    }
    if(is.null(group_vector)) {
      group_vector <- subset_group_vector
    } else {
      group_vector <- c(group_vector, subset_group_vector)
    }
    if(is.null(indiv_vector)) {
      indiv_vector <- subset_indiv_vector
    } else {
      indiv_vector <- c(indiv_vector, subset_indiv_vector)
    }
    if(is.null(plate_vector)) {
      plate_vector <- subset_plate_vector
    } else {
      plate_vector <- c(plate_vector, subset_plate_vector)
    }
    if(is.null(conc_vector)) {
      conc_vector <- subset_conc_vector
    } else {
      conc_vector <- c(conc_vector, subset_conc_vector)
    }
  }
  
  ns_all <- length(season_vector)

  season_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        season_kernel[i,j] <- 1
      } else if(season_vector[i] == "W" && season_vector[j] == "W") {
        season_kernel[i,j] <- within_group_corr
      } else if(season_vector[i] == "D" && season_vector[j] == "D") {
        season_kernel[i,j] <- within_group_corr
      } else {
        season_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  age_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        age_kernel[i,j] <- 1
      } else if(age_vector[i] == "J" && age_vector[j] == "J") {
        age_kernel[i,j] <- within_group_corr
      } else if(age_vector[i] == "A" && age_vector[j] == "A") {
        age_kernel[i,j] <- within_group_corr
      } else if(age_vector[i] == "G" && age_vector[j] == "G") {
        age_kernel[i,j] <- within_group_corr
      } else {
        age_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  group_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        group_kernel[i,j] <- 1
      } else if(group_vector[i] == group_vector[j]) {
        group_kernel[i,j] <- within_group_corr
      } else {
        group_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  indiv_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        indiv_kernel[i,j] <- 1
      } else if(indiv_vector[i] == indiv_vector[j]) {
        indiv_kernel[i,j] <- within_group_corr
      } else {
        indiv_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  plate_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        plate_kernel[i,j] <- 1
      } else if(plate_vector[i] == plate_vector[j]) {
        plate_kernel[i,j] <- within_group_corr
      } else {
        plate_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  conc_kernel <- matrix(0, nrow=ns_all, ncol=ns_all)
  for(i in 1:ns_all) {
    for(j in 1:ns_all) {
      if(i == j) {
        conc_kernel[i,j] <- 1
      } else if(conc_vector[i] == conc_vector[j]) {
        conc_kernel[i,j] <- within_group_corr
      } else {
        conc_kernel[i,j] <- between_group_corr
      }
    }
  }
  
  # residual (white noise kernel)
  if(include_residual) {
    empty_kernel <- matrix(within_group_corr, ns_all, ns_all)
    diag(empty_kernel) <- 1
  } else {
    empty_kernel <- matrix(0, ns_all, ns_all)
  }
  
  if(length(baboons) <= 0) {
    if("phyloseq" %in% class(data)) {
      tag <- ""
    } else {
      tag <- "metagenomics"
    }
    plot_cov(week_kernel, paste("date_sample_correlation", tag, sep=""))
    plot_cov(season_kernel, paste("season_sample_correlation", tag, sep=""))
    plot_cov(age_kernel, paste("age_sample_correlation", tag, sep=""))
    plot_cov(group_kernel, paste("group_sample_correlation", tag, sep=""))
    plot_cov(indiv_kernel, paste("indiv_sample_correlation", tag, sep=""))
    plot_cov(plate_kernel, paste("plate_sample_correlation", tag, sep=""))
    plot_cov(conc_kernel, paste("DNAconc_sample_correlation", tag, sep=""))
    if(include_residual) {
      plot_cov(empty_kernel, paste("residual_sample_correlation", tag, sep=""))
    }
  }

  # optimize the scale of each variance component
  logd_mat_t <- function(s) {
    return(logd_matrixt(s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8],
                        ilr_data, week_kernel, season_kernel, group_kernel, age_kernel,
                        indiv_kernel, plate_kernel, conc_kernel, empty_kernel))
  }
  
  components <- 8
  estimates <- matrix(0, components, optim_it)
  for(i in 1:optim_it) {
    cat("Optimization iteration",i,"\n")
    # about 3x faster if we have R optim call a compiled C++ function
    # could probably speed up further using DEoptim instead of optim here
    res <- optim(par=runif(components), fn=logd_mat_t, method="L-BFGS-B",
                 lower=rep(0.00001,components), upper=rep(10,components))
    estimates[,i] <- res$par
  }
  
  kernels <- c("weekly","seasonal","group","age","individual","plate","DNAconc","residual")
  if(optim_it > 1) {
    if(include_residual) {
      total_var <- colSums(estimates)
    } else {
      total_var <- colSums(estimates[1:7,])
    }
  } else {
    if(include_residual) {
      total_var <- c(sum(estimates))
    } else {
      total_var <- c(sum(estimates[1:7,]))
    }
  }
  for(i in 1:optim_it) {
    rank <- order(estimates[,i], decreasing=TRUE)
    cat("Iteration",i,":")
    for(j in 1:components) {
      if(rank[j] != components || include_residual) {
        cat(" ",kernels[rank[j]]," (",round(estimates[rank[j],i]/total_var[i],digits=3),")",sep="")
      }
    }
    cat("\n")
  }
}

# ====================================================================================================================
# DYNAMIC LINEAR MODELS
# ====================================================================================================================

# quick plot function for log ratios (x=time)
plot_lr <- function(ys, filename=NULL) {
  df <- gather_array(ys, "value", "time", "component")
  df$component <- as.factor(df$component)
  
  p <- ggplot(df, aes(x=time, y=value, group=component)) +
    geom_line(aes(color=component)) +
    theme_minimal() +
    theme(legend.position="none")
  if(!is.null(filename)) {
    ggsave(filename, scale=1.5, height=2, width=4)
  } else {
    show(p)
  }
}

# quick plot function for proportions (x=time)
plot_prop <- function(ys, filename=NULL) {
  # visualize as proportions
  proportions <- clrInv(ys)
  df <- gather_array(proportions, "value", "time", "component")
  df$component <- factor(df$component)
  
  categories <- unique(df$component)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df$component)))
  p <- ggplot(df, aes(x=time, y=value, fill=component)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none")
  if(!is.null(filename)) {
    ggsave(filename, scale=1.5, height=2, width=4)
  } else {
    show(p)
  }
}

# quick plot function for true thetas, filtered thetas, and smoothed thetas - sanity check
# data_obj is the output of five_taxa_simuation()
# fit_obj.f is the output of fit_filter()
# fit_obj.s is the output of fit_smoother()
plot_theta_fits <- function(data_obj, fit_obj.f=NULL, fit_obj.s=NULL, observation_vec=NULL, filename=NULL) {
  if(is.null(fit_obj.f) && is.null(fit_obj.s)) {
    return()
  }
  T <- nrow(data_obj$ys)
  T_actual <- T
  if(!is.null(observation_vec)) {
    T_actual <- max(observation_vec)
  }
  xs <- 1:T_actual
  D <- ncol(data_obj$ys)
  
  df.true <- gather_array(data_obj$ys, "logratio", "timepoint", "taxon")
  # replace timepoint with the actual day offset
  if(!is.null(observation_vec)) {
    for(i in 1:nrow(df.true)) {
      df.true$timepoint[i] <- observation_vec[df.true$timepoint[i]]
    }
  }
  df.true <- cbind(df.true, which="true")
  
  filtered.observations <- matrix(0, T_actual, D)
  smoothed.observations <- matrix(0, T_actual, D)
  for(t in 1:T_actual) {
    if(!is.null(fit_obj.f)) {
      filtered.observations[t,] <- data_obj$F%*%fit_obj.f$Thetas.t[,,t]
    }
    if(!is.null(fit_obj.s)) {
      smoothed.observations[t,] <- data_obj$F%*%fit_obj.s$Thetas.t[,,t]
    }
  }
  df.all <- df.true
  if(!is.null(fit_obj.f)) {
    df.filtered <- gather_array(filtered.observations, "logratio", "timepoint", "taxon")
    df.filtered <- cbind(df.filtered, which="filtered")
    df.all <- rbind(df.all, df.filtered)
  }
  if(!is.null(fit_obj.s)) {
    df.smoothed <- gather_array(smoothed.observations, "logratio", "timepoint", "taxon")
    df.smoothed <- cbind(df.smoothed, which="smoothed")
    df.all <- rbind(df.all, df.smoothed)
  }
  
  p <- ggplot(data=df.all, aes(x=timepoint, y=logratio, color=which, group=which)) + 
    geom_line() + 
    facet_wrap(~taxon) +
    theme_minimal()
  if(!is.null(filename)) {
    ggsave(filename, plot=p, scale=1.5, width=8, height=5)
  } else {
    show(p)
  }
}

# unfinished; busted
# fit_and_plot_intervals <- function(data_obj, observation_vec=observation_vec) {
#   S <- 10
#   T <- nrow(data_obj$ys)
#   if(!is.null(observation_vec)) {
#     T <- max(observation_vec)
#   }
#   draws <- array(dim=c(T, ncol(data_obj$ys), S)) # time (rows) x taxa (columns) x samples
#   for(s in 1:S) {
#     cat("Generating sample",s,"\n")
#     fit.f <- fit_filter(data_obj, observation_vec=observation_vec)
#     fit.s <- fit_smoother(data_obj, fit.f)
#     draws[,,s] <- t(fit.s$Thetas.t[1,,])
#   }
#   labeled_ys <- data_obj$ys
#   if(!is.null(observation_vec)) {
#     rownames(labeled_ys) <- observation_vec
#   }
#   samples <- gather_array(data_obj$ys, "coord", "timepoint", "taxon")
#   # crazy inefficient; fix this
#   samples <- cbind(samples, lower=NA)
#   samples <- cbind(samples, upper=NA)
#   for(i in 1:nrow(samples)) {
#     samples[i,"lower"] <- min(draws[samples[i,"timepoint"],samples[i,"taxon"],])
#     samples[i,"upper"] <- max(draws[samples[i,"timepoint"],samples[i,"taxon"],])
#   }
#   taxon <- 2
#   p <- ggplot(samples[samples$taxon==taxon,], aes(timepoint)) +
#     geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey80") +
#     geom_line(aes(y=coord))
#   p
# }

plot_mean_Sigma <- function(fit_obj, filename=NULL) {
  mean.Sigma <- fit_obj$Xi/(fit_obj$upsilon - nrow(fit_obj$Xi) - 1)
  if(!is.null(filename)) {
    png(filename)
    image(mean.Sigma)
    dev.off()
  } else {
    image(mean.Sigma)
  }
}


# generate Fourier form rotation matrix for a given omega (a function of period)
build_G <- function(period, harmonics=1) {
  if(harmonics == 1) {
    omega <- 2*pi/period
    return(matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2))
  } else {
    G <- matrix(0, 2*harmonics, 2*harmonics)
    for(h in 1:harmonics) {
      omega <- 2*pi*h/period
      common_offset <- (h-1)*2
      G[(common_offset+1),(common_offset+1)] <- cos(omega)
      G[(common_offset+1),(common_offset+2)] <- sin(omega)
      G[(common_offset+2),(common_offset+1)] <- -sin(omega)
      G[(common_offset+2),(common_offset+2)] <- cos(omega)
    }
    return(G)
  }
}

# 2-taxa simulation just for illustrative purposes
two_taxa_simulation <- function(T=30) {
  state_noise_scale <- 0.1
  gamma.t <- 1
  observational_noise_scale <- state_noise_scale*gamma.t
  
  period <- 10

  F <- matrix(c(1, 0), 1, 2)
  W.t <- diag(2)*state_noise_scale
  G <- build_G(period)
  upsilon <- 100
  
  # we'll do all combinations of noise and generate plots to see how the resultant log ratios
  # and proportions look
  noise_flags <- matrix(FALSE, 4, 2)
  noise_flags[2,1] <- TRUE # use state transition noise
  noise_flags[3,2] <- TRUE # use observational noise
  noise_flags[4,] <- TRUE # use both types of noise
  
  D <- 2
  
  for(i in 1:dim(noise_flags)[1]) {
    Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
    Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
    
    theta.t <- rmatrixnormal(1, matrix(1, 2, D), W.t, Sigma)[,,1]

    ys <- matrix(0, T, D)
    for(t in 1:T) {
      # state equation
      theta.t <- G%*%theta.t
      if(noise_flags[i,1]) {
        theta.t = theta.t + rmatrixnormal(1, matrix(0, 2, D), W.t, Sigma)[,,1]
      }
      # observation equation
      ys[t,] <- F%*%theta.t
      if(noise_flags[i,2]) {
        ys[t,] <- ys[t,] + rmvnorm(1, rep(0, D), gamma.t*Sigma)
      }
    }
    
    plot_lr(ys, filename=paste(i,"_1.png",sep=""))
    plot_prop(ys, filename=paste(i,"_2.png",sep=""))
  }
}

five_taxa_simulation <- function(T=40, indep_taxa=TRUE, uniform_start=TRUE, noise_scale=0.05, save_images=TRUE) {
  state_noise_scale <- noise_scale
  gamma <- 1
  observational_noise_scale <- state_noise_scale*gamma
  
  period <- 10

  F <- matrix(c(1, 0), 1, 2)
  W <- diag(2)*state_noise_scale
  G <- build_G(period)
  upsilon <- 100
  
  D <- 5
  
  if(indep_taxa) {
    # fully independent taxa
    Xi <- diag(D)*observational_noise_scale*(upsilon-D-1)
    tag <- "indeptax_"
  } else {
    # correlated/anti-correlated taxa
    Xi <- diag(D)
    Xi[1,2] <- 0.5
    Xi[2,1] <- Xi[1,2]
    Xi[3,4] <- -0.5
    Xi[4,3] <- Xi[3,4]
    Xi <- Xi*observational_noise_scale*(upsilon-D-1)
    tag <- "nonindeptax_"
  }
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  
  # very similar initial states; everybody is pretty much in phase here
  if(uniform_start) {
    M.0 <- matrix(1, 2, D)
    tag <- paste0(tag, "unifinit")
  } else {
    M.0 <- matrix(rnorm(10, 1, 1), 2, D)
    tag <- paste0(tag, "nonunifinit")
  }
  C.0 <- W
  
  theta.t <- rmatrixnormal(1, M.0, C.0, Sigma)[,,1]
  
  ys <- matrix(0, T, D)
  for(t in 1:T) {
    # state equation
    theta.t <- G%*%theta.t + rmatrixnormal(1, matrix(0, 2, D), W, Sigma)[,,1]
    # observation equation
    ys[t,] <- F%*%theta.t + rmvnorm(1, rep(0, D), gamma*Sigma)
  }
  
  if(save_images) {
    plot_lr(ys, filename=paste0("plots/DLMsim_logratios_",tag,".png"))
    plot_prop(ys, filename=paste0("plots/DLMsim_proportions_",tag,".png"))
  }
  
  return(list(ys=ys, F=F, W=W, G=G, upsilon=upsilon, Xi=Xi, gamma=gamma, Sigma=Sigma, M.0=M.0, C.0=C.0))
}

# powers G through eigenvalue decomposition
stack_G <- function(G, it_begin, it_end, descending=TRUE, transpose=FALSE) {
  obj <- G
  if(transpose) {
    obj <- t(G)
  }
  obj <- eigen(obj)
  e_vec <- obj$vectors
  e_val <- diag(2)
  diag(e_val) <- obj$values
  if(it_begin != it_end) {
    if(descending) {
      if(it_begin > it_end) {
        power_it <- length(1:(it_begin-it_end+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    } else {
      if(it_begin < it_end) {
        power_it <- length(1:(it_end-it_begin+1))
        e_val <- e_val**power_it
      } else {
        # invalid case
        e_val <- matrix(0, 2, 2)
      }
    }
  }
  # explicitly only returning the real part of A
  # some tiny complex eigenvalues can be produced -- cool to truncate in this way?
  ret_val <- Re(e_vec%*%e_val%*%solve(e_vec))
  return(ret_val)
}

build_A <- function(T, G, C.0, W, gamma, save_images=T) {
  # calculate A (covariance matrix over states) for this simulation
  # this is the expression exactly as in the manuscript, calculated from Cov(eta_t, eta_{t-k})
  A <- matrix(0, T, T)
  for(i in 1:T) {
    for(j in 1:T) {
      if(i == j) {
        # diagonal
        t <- j
        first_sum <- matrix(0, 2, 2)
        if(t >= 2) {
          for(ell in t:2) {
            G_left <- stack_G(G, t, ell)
            G_right <- stack_G(G, ell, t, descending=FALSE, transpose=TRUE)
            addend <- G_left%*%W%*%G_right
            first_sum <- first_sum + addend
          }
        }
        # second sum
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, t, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C.0%*%G_right
        A[t,t] <- gamma + F%*%(W + first_sum + second_sum)%*%t(F)
      } else {
        tk <- i
        t <- j
        if(j < i) {
          tk <- j
          t <- i
        }
        # off-diagonal
        first_sum <- matrix(0, 2, 2)
        for(ell in tk:2) {
          G_left <- stack_G(G, t, ell)
          G_right <- stack_G(G, ell, tk, descending=FALSE, transpose=TRUE)
          first_sum <- first_sum + G_left%*%W%*%G_right
        }
        G_left <- stack_G(G, t, 1)
        G_right <- stack_G(G, 1, tk, descending=FALSE, transpose=TRUE)
        second_sum <- G_left%*%C.0%*%G_right
        G_left <- stack_G(G, t, tk+1)
        A[i,j] <- F%*%(G_left%*%W + first_sum + second_sum)%*%t(F)
      }
    }
  }
  
  if(save_images) {
    png("plots/DLMsim_Amat.png")
    image(A)
    dev.off()
  }
  
  if(min(eigen(A)$values) < 0) {
    cat("Matrix A has negative eigenvalue(s)!\n")
  }
}

# Kalman filter
# data_obj is a list containing ys, F, W, G, upsilon, Xi, gamma, Sigma, M.0, C.0
# observation_vec (if present) indicates the spacing of observations, e.g. c(1, 3, 4, 7)
#   indicates observations 2, 5, & 6 are missing and should be imputed in the usual way
fit_filter <- function(data_obj, censor_vec=NULL, observation_vec=NULL, discount=NULL) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  if(!is.null(observation_vec)) {
    T <- max(observation_vec)
  } else {
    T <- nrow(data_obj$ys)
  }
  if(is.null(censor_vec)) {
    censor_vec <- rep(0, T)
  }
  upsilon.t <- data_obj$upsilon
  Xi.t <- data_obj$Xi
  M.t <- data_obj$M.0
  C.t <- data_obj$C.0
  Thetas.t <- array(0, dim=c(Theta.dim, D, T)) # sample at each t
  Cs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  Rs.t <- array(0, dim=c(Theta.dim, Theta.dim, T))
  W.t <- NULL # last observed W.t (when using discount)
  for(t in 1:T) {
    if(censor_vec[t] == 1 || (!is.null(observation_vec) && !(t %in% observation_vec))) {
      # case 1: no observations at this time point for any individuals
      # impute all by passing prior on as posterior
      P.t <- data_obj$G%*%C.t%*%t(data_obj$G)
      if(is.null(discount)) {
        R.t <- P.t + data_obj$W
      } else {
        # note: West & Harrison (p. 352) suggest not applying discount for missing observations
        if(is.null(W.t)) {
          W.t <- ((1-discount)/discount)*P.t
        }
        R.t <- P.t + W.t
      }
      R.t <- round(R.t, 10)
      Rs.t[,,t] <- R.t
      A.t <- data_obj$G%*%M.t
      # carry forward without update
      M.t <- A.t
      Ms.t[,,t] <- M.t
      C.t <- R.t
      Cs.t[,,t] <- C.t
      # no change to Sigma.t parameters
      Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
      Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
    } else {
      # case 2: apply update where possible, note: F.t.T is F for us here and G.t is G
      # one-step ahead forecast at t
      P.t <- data_obj$G%*%C.t%*%t(data_obj$G)
      if(is.null(discount)) {
        R.t <- P.t + data_obj$W
      } else {
        W.t <- ((1-discount)/discount)*P.t
        R.t <- P.t + W.t
      }
      R.t <- round(R.t, 10)
      Rs.t[,,t] <- R.t
      A.t <- data_obj$G%*%M.t
      f.t.T <- data_obj$F%*%A.t
      q.t <- data_obj$gamma + (data_obj$F%*%R.t%*%t(data_obj$F))[1,1]
      if(is.null(observation_vec)) {
        # we're assuming 1 individual if observation_mat is missing
        e.t.T <- data_obj$ys[t,] - f.t.T
      } else {
        # if more than one sample on the same day (mislabeled?) take the first for now
        e.t.T <- data_obj$ys[as(which(observation_vec == t)[1], "numeric"),] - f.t.T
      }
      S.t <- R.t%*%t(data_obj$F)/q.t
      M.t <- A.t + S.t%*%e.t.T
      Ms.t[,,t] <- M.t
      C.t <- R.t - q.t*S.t%*%t(S.t)
      C.t <- round(C.t, 10)
      Cs.t[,,t] <- C.t
      upsilon.t <- upsilon.t + 1
      Xi.t <- Xi.t + t(e.t.T)%*%e.t.T/q.t
      Sigma.t <- rinvwishart(1, upsilon.t, Xi.t)[,,1]
      Thetas.t[,,t] <- rmatrixnormal(1, M.t, C.t, Sigma.t)[,,1]
    }
  }
  return(list(Thetas.t=Thetas.t, upsilon=upsilon.t, Xi=Xi.t, Ms.t=Ms.t, Cs.t=Cs.t, Rs.t=Rs.t))
}

# simulation smoother
# fit_obj is the output of fit_filter()
fit_smoother <- function(data_obj, fit_obj) {
  D <- ncol(data_obj$ys)
  Theta.dim <- ncol(data_obj$G)
  T <- dim(fit_obj$Thetas.t)[3]
  Sigma.t <- rinvwishart(1, fit_obj$upsilon, fit_obj$Xi)[,,1]
  Thetas.t.smoothed <- array(0, dim=c(Theta.dim, D, T))
  Ms.t <- array(0, dim=c(Theta.dim, D, T))
  etas.t <- matrix(0, D, T)
  rmatnorm_mean <- fit_obj$Ms.t[,,T]
  dim(rmatnorm_mean) <- c(Theta.dim, D) # fix if the state has one dimension only
                                        # can't do drop=FALSE here
  Thetas.t.smoothed[,,T] <- rmatrixnormal(1, rmatnorm_mean, fit_obj$Cs.t[,,T], Sigma.t)[,,1]
  etas.t[,T] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,T], 1, Sigma.t)[,,1]
  Ms.t[,,T] <- rmatnorm_mean
  for(t in (T-1):1) {
    Z.t <- fit_obj$Cs.t[,,t]%*%t(data_obj$G)%*%solve(fit_obj$Rs.t[,,(t+1)])
    M.t.star <- fit_obj$Ms.t[,,t] + Z.t%*%(Thetas.t.smoothed[,,(t+1)] - data_obj$G%*%fit_obj$Ms.t[,,t])
    Ms.t[,,t] <- M.t.star
    C.t.star <- round(fit_obj$Cs.t[,,t] - Z.t%*%fit_obj$Rs.t[,,(t+1)]%*%t(Z.t), 10)
    Thetas.t.smoothed[,,t] <- rmatrixnormal(1, M.t.star, C.t.star, Sigma.t)[,,1]
    # draw an eta, this is how we'll estimate modeled covariance
    etas.t[,t] <- rmatrixnormal(1, data_obj$F%*%Thetas.t.smoothed[,,t], 1, Sigma.t)[,,1]
  }
  return(list(Thetas.t=Thetas.t.smoothed, etas.t=etas.t, Ms.t=Ms.t))
}

pull_indiv_data <- function(sname, asv_data, pseudocount=0.5, date_begin="1900-01-01", date_end="2100-01-01",
                               subset_dim=0, lr_transform=TRUE) {
  pruned <- prune_samples(sample_data(non_reps)$sname==sname, non_reps)
  pruned <- prune_samples((sample_data(pruned)$collection_date > date_begin) &
                            (sample_data(pruned)$collection_date < date_end), pruned)
  cat(paste0("Real data set (",sname,") has ",nsamples(pruned)," samples and ",ntaxa(pruned)," taxa\n"))
  counts <- otu_table(pruned)@.Data # samples (rows) x taxa (columns)
  if(lr_transform) {
    ys <- driver::clr(counts + pseudocount)
  } else {
    ys <- counts
  }
  # the baseline date common to individuals
  min_date <- min(sample_data(pruned)$collection_date)
  dates_observed <- sample_data(pruned)$collection_date
  observation_vec <- sapply(dates_observed, function(x) { round(difftime(x, min_date, units="days"))+1 } )
  # subset taxa
  if(subset_dim > 0) {
    ys <- ys[,1:subset_dim]
  }
  return(list(ys=ys, observation_vec=observation_vec, sname=sname))
}

# ====================================================================================================================
# GAUSSIAN PROCESS WORKUP : KERNELS
# ====================================================================================================================

# X is Q x N as in other kernels
# bandwidth: rho as chosen gives antiphase observations a correlation of ~0.1
PER <- function(X, sigma=1, rho=1, period=24, jitter=1e-10){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

WHITENOISE <- function(X, sigma=1, jitter=1e-10) {
  dist <- as.matrix(dist(t(X)))
  G <- diag(ncol(dist))*sigma^2 + jitter*diag(ncol(dist))
  return(G)
}











