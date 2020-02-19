# this file collects sample IDs from saved PAM (partitioning around medoids) output; these are prospective
# metagenomics samples with near-optimal spread in terms of their variance

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

md <- readRDS(file.path(relative_path,output_dir,"PAM_filtered_metadata.rds"))
for(k in c(500, 1000, 1500, 2000)) {
  load(file.path(relative_path,output_dir,paste0("PAM_clustering_",k,".RData")))
  df <- data.frame(SID=md$sample_id[sample_cr$medoids],
                   collection_date=md$collection_date[sample_cr$medoids],
                   sname=md$sname[sample_cr$medoids],
                   sample_status=md$sample_status[sample_cr$medoids])
  write.table(df, file=file.path(relative_path,output_dir,paste0("PAM_selected_samples_",k,".tsv")), quote=FALSE, sep='\t', col.names = NA)
}

