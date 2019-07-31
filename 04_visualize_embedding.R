library(Rcpp)
library(ggplot2)
library(dplyr)
source("include.R")

# allowable label types: individual, group, counts, density

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) {
  stop("Testing usage: Rscript 04_visualize_embedding.R family Lambda counts", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
label_type <- args[3]

get_other_labels <- function(df, data, individuals, annotation="group") {
  labels <- numeric(nrow(df))
  names(labels) <- df$label
  if(annotation == "group") {
    metadata <- sample_data(data)
    primary_group <- metadata %>%
        select(c("sname", "collection_date", "grp")) %>%
        filter(sname %in% individuals) %>% 
        group_by(sname, grp) %>%
        tally() %>%
        slice(which.max(n))
    for(indiv in individuals) {
      labels[df$label == indiv] <- primary_group[primary_group$sname == indiv,]$grp[[1]]
    }
  }
  if(annotation == "matgroup") {
    metadata <- sample_data(data)
    for(indiv in individuals) {
      labels[df$label == indiv] <- metadata[metadata$sname == indiv,]$matgrp[[1]]
    }
  }
  if(annotation == "counts" | annotation == "density") {
    for(indiv in individuals) {
      cat("Parsing individual",indiv,"...\n")
      indiv <<- indiv # needs to be global (bug)
      indiv_subset <- subset_samples(data, sname==indiv)
      sample_count <- phyloseq::nsamples(indiv_subset)
      labels[df$label == indiv] <- round(sample_count, -1) # discretize
    }
  }
  if(annotation == "density") {
    metadata <- sample_data(data)
    md_subset <- metadata %>%
        select(c("sname", "collection_date")) %>%
        filter(sname %in% individuals) %>%
        group_by(sname) %>%
        summarize(delta=difftime(max(collection_date), min(collection_date), units="days"))
    for(indiv in individuals) {
      labels[names(labels) == indiv] <- labels[names(labels) == indiv]/md_subset[md_subset$sname == indiv,]$delta[[1]]
    }
    labels <- round(labels*100)
  }
  df2 <- data.frame(x=df$x, y=df$y, z=df$z, labels=as.factor(labels))
  return(df2)
}

plot_axes <- function(df, df_centroids=NULL, axis1="x", axis2="y", label_type="individual", legend=TRUE) {
  p <- ggplot() + geom_point(data=df, aes_string(x=axis1, y=axis2, color="labels"))
  if(label_type == "individual") {
    # label the centroids directly
    p <- p + geom_text(data=df_centroids, aes_string(x=paste0("mean_",axis1), y=paste0("mean_",axis2), label="labels"), color="black", fontface="bold")
  }
  if(!legend | label_type == "individual") {
    p <- p + theme(legend.position='none')
  }
  plot_save_name <- paste0(which_measure,"_ordination_",label_type,"_",axis1,axis2,".png")
  img_width <- 4
  if(legend) {
    img_width <- 4.5
  }
  ggsave(paste0("plots/basset/",level,"/",plot_save_name), plot=p, scale=2,
           width=img_width, height=4, units="in", dpi=100)
}

plot_ordination <- function(level, which_measure, label_type, legend=TRUE) {
  load(paste0("plots/basset/",level,"/",which_measure,"_ordination.RData")) # 'df' object
  load(paste0("plots/basset/",level,"/",which_measure,"_ordination_centroids.RData")) # 'df_centroids' object
  if(label_type != "individual") {
    glom_data <- load_glommed_data(level=level, replicates=TRUE)
    df <- get_other_labels(df, glom_data, unique(df$label), annotation=label_type)
  }
  # axes 1 & 2
  plot_axes(df, df_centroids, "x", "y", label_type, legend=legend)
  # axes 2 & 3
  plot_axes(df, df_centroids, "y", "z", label_type, legend=legend)
  # axes 1 & 3
  plot_axes(df, df_centroids, "x", "z", label_type, legend=legend)
}

plot_ordination(level, which_measure, label_type)
