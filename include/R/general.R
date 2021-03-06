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

data_dir <- "input"
output_dir <- "output"
model_dir <- "output/model_fits"
plot_dir <- "output/plots"

pc <- 0.5

base_path <- "/data/mukherjeelab/rulesoflife"
#base_path <- "C:/Users/kim/Documents/rules_of_life"

sourceCpp(file.path(base_path,"include/cpp/fast_corr.cpp"))

# these lists were generated manually
# 10 max-sampled individuals
best_sampled <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI")
# individuals with # samples >= 100
over_100 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
              "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO")
over_75 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
             "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO", "LAD", "LOB", "WAD", "GAB", "LIW",
             "VIN", "TAL", "VEX", "VEI", "ALE", "MBE", "WHE", "WYN", "LOL", "HOL", "NOB", "VOT", "LYE",
             "HON", "DAG", "DUN", "OTI", "LUI", "OFR", "LAZ", "ONY", "VEL", "ELV", "FAX", "ORI", "EAG",
             "ODE", "NIK", "VAP", "WIP")
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
# individuals with # samples >= 40
over_40 <- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
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
             "ETO", "KEL", "NUT", "WES", "IDI", "ISO", "PRU", "YOG", "ZAI")
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

# designate which of the individual/sname subsets to use throughout!
sname_list <- over_40

# ====================================================================================================================
# ACCESSORS
# ====================================================================================================================

# read unagglomerated, unfiltered 16S data
read_data <- function(replicates=TRUE) {
  # rows are samples, columns are taxa
  if(replicates) {
    cat("Using replicates...\n")
    data <- readRDS(file.path(base_path,data_dir,"original_16S/emp_baboon_pool_T_w_techReps.RDS"))
  } else {
    cat("NOT using replicates...\n")
    data <- readRDS(file.path(base_path,data_dir,"original_16S/emp_baboon_NewFiltr.RDS"))
  }
  # for now, just remove non-Bacterial domain
  data <- subset_taxa(data, domain=="Bacteria")
  return(data)
}

# pull metadata from phyloseq object
read_metadata <- function(data) {
  # metadata: rows are samples, "sample-specific variables" are columns
  metadata <- phyloseq::sample_data(data)
  # sort ASC by collection date
  metadata <- metadata[order(metadata$sname,metadata$collection_date),]
  return(metadata)
}

# load agglomerated data
load_glommed_data <- function(level="species", replicates=TRUE) {
  if(replicates) {
    filename <- file.path(base_path,data_dir,paste0("glom_data_",level,"_reps_tree.rds"))
    if(file.exists(filename)) {
      glom_data <- readRDS(filename)
    } else {
      stop(paste0("Agglomerated data file for ",level," (+replicates) does not exist:",filename))
    }
  } else {
    filename <- file.path(base_path,data_dir,paste0("glom_data_",level,".rds"))
    if(file.exists(filename)) {
      glom_data <- readRDS(filename)
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

# filters `data` to exclude Mitochondrial and Chloroplast counts and puts low counts (below the thresholds
# specified by `count_threshold` and `sample_threshold`) into an "Other" category
filter_data <- function(data, level="genus", count_threshold=5, sample_threshold=20, verbose=TRUE) {
  cat("Filtering data at level",level,"with count threshold",count_threshold,"and sample threshold",sample_threshold,"\n")
  filter_fn <- file.path(base_path,data_dir,paste0("filtered_",level,"_",count_threshold,"_",sample_threshold,".rds"))
  if(!file.exists(filter_fn)) {
    snames <- unique(sample_data(data)$sname)
    counts <- otu_table(data)@.Data
    total_counts <- sum(counts)
    # a taxon must appear in at least {count_threshold} abundance in at least {sample_threshold} samples per individual
    keep_indices <- rep(TRUE, ntaxa(data))
    for(sname in snames) {
      cat("Parsing individual:",sname,"\n")
      sname <<- sname
      indiv_data <- subset_samples(data, sname==sname)
      # these are the indices to remove!
      keep_indices_indiv <- apply(counts, 2, function(x) sum(x >= count_threshold)/phyloseq::nsamples(indiv_data) >= sample_threshold/100)
      keep_indices <- keep_indices & keep_indices_indiv
    }
    collapse_indices <- !keep_indices
    saveRDS(collapse_indices, file.path(base_path,data_dir,paste0("filtered_idx_",level,"_",count_threshold,"_",sample_threshold,".rds")))
    # collapse mitochondria too
    tt <- tax_table(data)@.Data
    collapse_indices[which(tt[,colnames(tt) == "family"] == "Mitochondria")] <- TRUE
    collapse_indices[which(tt[,colnames(tt) == "order"] == "Chloroplast")] <- TRUE
    merged_data <- merge_taxa(data, which(collapse_indices == TRUE), 1)
    retained_counts <- sum(counts[,!collapse_indices])
    if(verbose) {
      cat("Collapsing",sum(collapse_indices == TRUE),"taxa of",ntaxa(data),"\n")
      collapse_tidx <- which(collapse_indices == TRUE)[1]
      cat("\tOther category collapsed into index:",collapse_tidx,sep=" ","\n")
      cat("\tTaxonomy of collapsed (BEFORE):",tax_table(merged_data)@.Data[collapse_tidx,],"\n")
      tax_table(merged_data)@.Data[collapse_tidx,] <- rep("Collapsed",7)
      cat("\tTaxonomy of collapsed (AFTER):",tax_table(merged_data)@.Data[collapse_tidx,],"\n")
      cat("\tCollapsed counts:",(total_counts-retained_counts),"of",total_counts,"(",(total_counts-retained_counts)/total_counts,"total )\n")
      cat("\tPercent zero-count in data set:",(1 - get_tiny_counts(merged_data, 1)),"\n")
    } else {
      tax_table(merged_data)@.Data[collapse_tidx,] <- rep("Collapsed",7)
    }
    saveRDS(merged_data, file=filter_fn)
  } else {
    merged_data <- readRDS(filter_fn)
  }
  return(merged_data)
}

load_and_filter <- function(level="family") {
  if(level == "ASV" | is.null(level)) {
    level <- "ASV"
    glom_data <- read_data(replicates=FALSE)
  } else {
    glom_data <- load_glommed_data(level=level, replicates=TRUE)
  }
  data <- filter_data(glom_data, level=level, verbose=FALSE)
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

# assumes samples x taxa
pick_alr_ref <- function(data) {
  counts <- otu_table(data)@.Data
  median <- which(order(apply(counts, 2, mean)) == round(ncol(counts)/2))
  return(median)
}

default_ALR_prior <- function(D, log_var_scale=1) {
  upsilon <- D-1+10 # lesser certainty
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*log_var_scale)%*%t(GG) # take diag as covariance over log abundances
  Xi <- Xi*(upsilon-D-1)
  return(list(upsilon=upsilon, Xi=Xi))
}
