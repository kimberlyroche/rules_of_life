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
  if(level == "species") {
    if(replicates) {
      if(file.exists("glom_data_species_reps.RData")) {
        load("glom_data_species_reps.RData")
      } else {
        stop("Agglomerated data file foe species (+replicates) does not exist.")
      }
    } else {
      if(file.exists("glom_data_species.RData")) {
        load("glom_data_species.RData")
      } else {
        stop("Agglomerated data file foe species (-replicates) does not exist.")
      }
    }
  }
  if(level == "genus") {
    if(replicates) {
      if(file.exists("glom_data_genus_reps.RData")) {
        load("glom_data_genus_reps.RData")
      } else {
        stop("Agglomerated data file foe genus (+replicates) does not exist.")
      }
    } else {
      if(file.exists("glom_data_genus.RData")) {
        load("glom_data_genus.RData")
      } else {
        stop("Agglomerated data file foe genus (-replicates) does not exist.")
      }
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
filter_data <- function(data, count_threshold=3, sample_threshold=0.2, verbose=FALSE) {
  total_counts <- sum(otu_table(data)@.Data)
  # get indices to collapse
  retained_counts <- sum(otu_table(filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= sample_threshold, prune=T))@.Data)
  collapse_indices <- !as.logical(filter_taxa(data, function(x) sum(x >= count_threshold)/phyloseq::nsamples(data) >= sample_threshold, prune=F))
  collapse_indices <- which(collapse_indices)
  merged_data <- merge_taxa(data, collapse_indices, 1)
  if(verbose) {
    cat("Collapsing",length(collapse_indices),"taxa of",ntaxa(data),"\n")
    cat("\tOther category collapsed into:",as.vector(tax_table(data)[collapse_indices[1]]),sep=" ","\n")
    cat("\tOther sequence variant:",rownames(tax_table(data))[collapse_indices[1]],sep=" ","\n")
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

# read in PiPhillin data from .tsv file
# requires 16S metadata so we can pare down the samples to those in-common between data sets
#
# returns a matrix where rows (and rownames) are enzymes and columns (and column names) are
# samples ordered by collection date (but not sub-ordered by anything)
read_metagenomics <- function(metadata) {
  piphillin_dir <- "original_data/Piphillin_20190222"
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
  idx.subset <- sname_subset_idx(metadata, sname)
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
  # df.prop.ordered$sample <- factor(df.prop.ordered$sample, levels=unique(df.prop.ordered$sample))
  # df.prop.ordered$enzyme <- factor(df.prop.ordered$enzyme, levels=unique(df.prop.ordered$enzyme))
  return(df.prop.ordered)
}

# expects a tidy array with columns |enzyme{factor;number}| |sample{factor;date}| |proportion{float}}
plot_timecourse_metagenomics <- function(metagenomics_prop, save_filename="metagenomics_timecourse", gapped=FALSE, legend=FALSE) {
  na.string <- ".N/A"

  # make sure the dates are in order and fix the order by converting to factors
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
  }
  df2 <- df2[order(df2$sample),]
  df2$enzyme <- factor(df2$enzyme, levels=unique(df2$enzyme))
  
  categories <- unique(df2$enzyme)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df2$enzyme)))
  if(gapped) {
    img_width <- 15
  } else {
    img_width <- 10
  }
  p <- ggplot(df2, aes(x=sample, y=proportion, fill=enzyme)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8))
  if(legend) {
    p <- p + theme(legend.position="bottom")
  } else {
    p <- p + theme(legend.position="none")
  }
  ggsave(paste("plots/", save_filename,".png",sep=""), plot=p, scale=2, width=img_width, height=4, units="in")
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
plot_timecourse_phyloseq <- function(data, save_filename, gapped=FALSE, legend=TRUE, legend_level="family") {
  n <- nsamples(data)
  p <- apply_proportion(data)
  df <- psmelt(p)
  df2 <- bind_cols(list(OTU=df$OTU, Sample=df$Sample, Abundance=df$Abundance))

  # for gapped plots, this is the OTU ID placeholder we'll use for dummy data
  na.string <- ".N/A"
  
  # replace Sample ID's with their dates for readability
  metadata <- sample_data(data)
  for(i in 1:dim(df2)[1]) {
    df2$Sample[i] <- metadata[metadata$sample_id==df2$Sample[i],"collection_date"][[1]]
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
  # replace ASV sequences with their (abbrev.) taxonomy for readability
  for(i in 1:dim(df3)[1]) {
    # show labels as order/family/genus
    # species is NA for all
    if(df3$OTU[i] != na.string) {
      if(legend_level == "species") {
        df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4:6]),collapse="/")
      } else if(legend_level == "genus") {
        df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4:5]),collapse="/")
      } else {
        df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4]),collapse="/")
      }
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
  p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8))
  if(n < 20) {
    # these are likely replicates
    img_width <- 5
  }
  if(legend) {
    p <- p + theme(legend.position="bottom")
  } else {
    p <- p + theme(legend.position="none")
  }
  ggsave(paste("plots/",save_filename,".png",sep=""), plot=p, scale=2, width=img_width, height=4, units="in")
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
                                 resample_rate=0.2) {
  if(date_diff_units != "weeks" && date_diff_units != "months" && date_diff_units != "seasons") {
    date_diff_units <- "weeks"
  }

  individuals <- unique(metadata$sname)
  individuals <- individuals[1:20] # for testing

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
  
  log_ratios <- apply_ilr(data)
  log_ratios <- scale(log_ratios, center=T, scale=F)
  for(r in 1:rounds) {
    if(resample) {
      cat("Resampling iteration",r,"\n")
    }
    lag.sums <- numeric(lag.max)
    lag.measured <- numeric(lag.max)
    for(indiv in individuals) {
      # pick an individual
      # this weird syntactic hack seems to be necessary for subset_samples?
      # apparently the thing you're filtering against must be globally available
      indiv <<- indiv
      # pull sample IDs associated with this individual
      sample_info <- metadata[metadata$sname == indiv, c("sample_id", "collection_date")]
      # appear already ordered but let's be paranoid
      sample_info <- sample_info[order(sample_info$collection_date),]
      # pull log ratios associated with this individual
      # there may be cases where none of the identified samples are present in this dataset
      # catch that case
      if(length(intersect(rownames(log_ratios), sample_info$sample_id)) > 0) {
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
          lags[lag,r] <- lag.sums[lag]/lag.measured[lag]
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

plot_mean_autocorrelation <- function(lags, filename="autocorrelation") {
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
  ggsave(paste(filename,".png",sep=""), plot=p, scale=2, width=4, height=3, units="in")
}

plot_bounded_autocorrelation <- function(lags, filename="autocorrelation") {
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
  ggsave(paste(filename,".png",sep=""), plot=p, scale=2, width=4, height=3, units="in")
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
  
  if(length(baboons) <= 25) {
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
  if(include_residual) {
    total_var <- colSums(estimates)
  } else {
    total_var <- colSums(estimates[1:7,])
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
