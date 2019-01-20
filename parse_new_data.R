# taxonomic information
tax <- readRDS("original_data/baboon_project_dada2_seqtab_preset_pooling_taxID.rds")
cat("Taxonomy file:",dim(tax)[1],"taxa with",dim(tax)[2],"level taxonomy info\n")

# sequence data (refiltered 1/16/19)
seq <- readRDS("original_data/merged_baboon_data_with_metadata_preset_more_blanks_greater_1000.rds")
cat("Sequence data file:",dim(seq)[1],"samples of",dim(seq)[2],"taxa\n")

# filter out (1) samples with all zeros (2) taxa will all zeros
keep_taxa <- colSums(seq)>0
keep_samples <- rowSums(seq)>0
seq <- seq[,keep_taxa]
seq <- seq[keep_samples,]
tax <- tax[keep_taxa,]
cat("Taxonomy file (all-zero filtered):",dim(tax)[1],"taxa with",dim(tax)[2],"level taxonomy info\n")
cat("Sequence data file (all-zero filtered):",dim(seq)[1],"samples of",dim(seq)[2],"taxa\n")
saveRDS(tax, "original_data/merged_baboon_tax_filtered.rds")
saveRDS(seq, "original_data/merged_baboon_seq_filtered.rds")

# sample metadata still in emp_baboon.rds; need to request new/full version?
phyloseq_obj <- readRDS("original_data/emp_baboon.RDS")
sample_metadata <- phyloseq::sample_data(phyloseq_obj)
cat("Sample metadata has",dim(sample_metadata)[1],"samples and",dim(sample_metadata)[2],"attributes\n")

# how many samples in seq have metadata present in the (old) metadata file?
sample_ids <- rownames(seq)
cat("There are",nsamples(sample_metadata[sample_metadata$sample_id %in% sample_ids,]),"of",
    dim(seq)[1],"samples with metadata present in the old metadata file\n")

# filter to these (temporarily)
shared_sample_ids <- as.vector(unlist(sample_metadata[sample_metadata$sample_id %in% sample_ids,"sample_id"]))
filtered_seq <- seq[shared_sample_ids,]
saveRDS(filtered_seq, "original_data/merged_baboon_seq_filtered_temp.rds")




# processed
tax <- readRDS("original_data/merged_baboon_tax_filtered.rds")
cat("Taxonomy file:",dim(tax)[1],"taxa with",dim(tax)[2],"level taxonomy info\n")

# sequence data (refiltered 1/16/19)
seq <- readRDS("original_data/merged_baboon_seq_filtered_temp.rds")
cat("Sequence data file:",dim(seq)[1],"samples of",dim(seq)[2],"taxa\n")
