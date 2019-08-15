library(phyloseq)
library(stray)
library(dplyr)
library(driver)
library(ggplot2)
library(scales)
library(coda)
library(Rcpp)
library(RcppEigen)
library(LaplacesDemon) # various sampling distributions
library(RColorBrewer)
library(tidyverse)
library(stringr)

data_dir <- "data/"
model_dir <- "output/model_fits/"
plot_dir <- "output/plots/"

pc <- 0.5

sourceCpp("include/cpp/fast_corr.cpp")

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
    filename <- paste0(data_dir,"glom_data_",level,"_reps.RData")
    if(file.exists(filename)) {
      load(filename)
    } else {
      stop(paste0("Agglomerated data file for ",level," (+replicates) does not exist:",filename))
    }
  } else {
    filename <- paste0(data_dir,"glom_data_",level,".RData")
    if(file.exists(filename)) {
      load(filename)
    } else {
      stop(paste0("Agglomerated data file for ",level," (-replicates) does not exist:",filename))
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

# filter data below a count threshold in a minimum number of samples into an "Other" category
# e.g. count_threshold=3, sample_threshold=0.2 filters taxa with no more than a 2-count in 80% of samples into an
# "Other" category, labeled by an arbitrary sequence variant in that "Other" category (see the print statement
# below identifying it for reference)
filter_data <- function(data, count_threshold=10, sample_threshold=0.66, collapse_level=NULL, verbose=TRUE) {
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

load_and_filter <- function(level="family") {
  glom_data <- load_glommed_data(level=level, replicates=TRUE)
  data <- filter_data(glom_data, collapse_level=level, verbose=FALSE)
  return(data)
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

# print human-readable taxonomy names; useful for figure labeling
print_readable_taxonomy <- function(level="family") {
  load_and_filter(level=level)
  tt <- tax_table(data)
  rownames(tt) <- NULL
  for(i in 1:(nrow(tt)-1)) {
    cat(tt[i,5],"/",tt[nrow(tt),5],"\n")
  }
}

# stray uses the D^th element as the ALR reference by default
# do some row shuffling in Y to put the reference at the end if there's
# some specific reference we want to use
reorient_count_matrix <- function(Y, alr_ref) {
  Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  return(Y)
}

default_ALR_prior <- function(D, log_var_scale=1) {
  upsilon <- D-1+10 # lesser certainty
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*log_var_scale)%*%t(GG) # take diag as covariance over log abundances
  Xi <- Xi*(upsilon-D-1)
  return(list(upsilon=upsilon, Xi=Xi))
}
