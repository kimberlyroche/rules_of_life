# this file embeds the (already calculated) posterior distances via UMAP

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

library(umap)

level <- "family"
n_samples <- 100
d_MAP <- readRDS(file.path(relative_path,output_dir,paste0("saved_distance_Sigma_",level,"_MAP.rds")))
d_full <- readRDS(file.path(relative_path,output_dir,paste0("saved_distance_Sigma_",level,".rds")))

data <- load_and_filter(level)
md <- sample_data(data)

# embed MAP
embedding_MAP <- umap(d_MAP, config=umap.defaults, method="naive")

# get group labels
group_labels <- list()
for(host in sname_list) {
  sub_samples <- subset_samples(data, sname == host)
  groups <- as.data.frame(sample_data(sub_samples))$grp
  grouped_groups <- table(groups)[1]
  primary_group <- names(which(grouped_groups == max(grouped_groups)))
  group_labels[host] <- primary_group
}

df <- data.frame(x=embedding$layout[,1], y=embedding$layout[,2], group=as.factor(unlist(group_labels)))
ggplot(df, aes(x=x, y=y, color=group)) +
  geom_point()

# embed full; this takes a few minutes
embedding_full <- umap(d_full, config=umap.defaults, method="naive")

host_labels <- c()
for(host in sname_list) {
  host_labels <- c(host_labels, rep(host, n_samples))
}

# this is crazy but confirms the main points that individuals are:
#    (1) fully distinct w/r/t their posteriors and 
#    (2) continuously different
df <- data.frame(x=embedding_full$layout[,1], y=embedding_full$layout[,2], host=as.factor(host_labels))
ggplot(df, aes(x=x, y=y, color=host)) +
  geom_point()
