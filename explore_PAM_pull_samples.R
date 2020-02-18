md <- readRDS("output/PAM_filtered_metadata.rds")
for(k in c(500, 1000, 1500, 2000)) {
  load(paste0("output/PAM_clustering_",k,".RData"))
  df <- data.frame(SID=md$sample_id[sample_cr$medoids],
                   collection_date=md$collection_date[sample_cr$medoids],
                   sname=md$sname[sample_cr$medoids],
                   sample_status=md$sample_status[sample_cr$medoids])
  write.table(df, file=paste0("output/PAM_selected_samples_",k,".tsv"), quote=FALSE, sep='\t', col.names = NA)
}

