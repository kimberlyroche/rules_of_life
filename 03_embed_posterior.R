library(Rcpp)
library(ggplot2)
library(dplyr)
source("include.R")

plot_Sigma_ordination <- function(fit, labels, label_type, legend=TRUE) {
  df <- data.frame(x=fit$points[,1], y=fit$points[,2], labels=labels)
  df_centroids <- df %>% group_by(labels) %>% summarise(mean_x=mean(x), mean_y=mean(y))
  save(df_centroids, file=paste0("plots/basset/",level,"/Sigma_centroids_",label_type,".RData"))
  p <- ggplot() +
         geom_point(data=df, aes(x=x, y=y, color=labels))
  if(!legend) {
    # label the centroids directly
    p <- p + geom_text(data=df_centroids, aes(x=mean_x, y=mean_y, label=labels), color="black", fontface="bold")
  }
  if(use_Riemann) {
    p <- p + ggtitle("Riemannian distance")
  } else {
    p <- p + ggtitle("Frobenius norm of difference")
  }
  if(!legend) {
    p <- p + theme(legend.position='none')
  }
  plot_save_name <- paste0("Sigma_posterior_ordination_",label_type,"_")
  if(use_Riemann) {
    plot_save_name <- paste0(plot_save_name,"Riemann.png")
  } else {
    plot_save_name <- paste0(plot_save_name,"Frobenius.png")
  }
  ggsave(paste0("plots/basset/",level,"/",plot_save_name), scale=2,
           width=4, height=4, units="in", dpi=100)
}

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Testing usage: Rscript 03_embed_posterior.R family 26", call.=FALSE)
}
level <- args[1]
D <- as.numeric(args[2])
use_Riemann <- TRUE

# testing
#date_lower_limit <- "2001-10-01"
#date_upper_limit <- "2003-11-30"
date_lower_limit <- NULL
date_upper_limit <- NULL

sourceCpp("cov_viz_test.cpp")

fitted_models <- list.files(path=paste0("subsetted_indiv_data/",level), pattern="*_bassetfit_.RData", full.names=TRUE, recursive=FALSE)
individuals <- sapply(fitted_models, function(x) { idx <- regexpr("_bassetfit_.RData", x); return(substr(x, idx-3, idx-1)) } )
names(individuals) <- NULL

# get group membership; use social group in which this individual spent the largest time
glom_data <- load_glommed_data(level=level, replicates=TRUE)
metadata <- read_metadata(glom_data)

if(!is.null(date_lower_limit) & !is.null(date_upper_limit)) {
  primary_group <- metadata %>%
    select(c("sname", "collection_date", "grp")) %>%
    filter(sname %in% individuals) %>% 
    filter(collection_date >= date_lower_limit) %>%
    filter(collection_date <= date_upper_limit) %>%
    group_by(sname, grp) %>%
    tally() %>%
    slice(which.max(n))
} else {
  primary_group <- metadata %>%
    select(c("sname", "collection_date", "grp")) %>%
    filter(sname %in% individuals) %>% 
    group_by(sname, grp) %>%
    tally() %>%
    slice(which.max(n))
}

cat("Group membership:\n")
print(tbl_df(primary_group), n=50)

n_samples <- 100 # eta samples to draw; Kalman filter, simulation smoother run on these

n_samples_subset <- 100

n_indiv <- length(individuals)
all_samples <- matrix(NA, D-1, (D-1)*n_samples_subset*n_indiv)
indiv_labels <- c()
group_labels <- c()
for(i in 1:n_indiv) {
  load(paste0("subsetted_indiv_data/",level,"/",individuals[i],"_bassetfit_.RData"))
  Sigma <- Sigma[,,1:n_samples_subset]
  all_samples[,((i-1)*(D-1)*n_samples_subset+1):(i*(D-1)*n_samples_subset)] <- Sigma
  indiv_labels <- c(indiv_labels, rep(individuals[i], n_samples_subset))
  group_labels <- c(group_labels, rep(primary_group[primary_group$sname == individuals[i], c("grp")][[1]], n_samples_subset))
}
group_labels <- as.factor(group_labels)

distance_mat <- matrix(NA, n_samples_subset*n_indiv, n_samples_subset*n_indiv)
for(i in 1:(n_indiv*n_samples_subset)) {
  for(j in i:(n_indiv*n_samples_subset)) {
    i_idx <- (i-1)*(D-1)
    A <- all_samples[,(i_idx+1):(i_idx+(D-1))]
    j_idx <- (j-1)*(D-1)
    B <- all_samples[,(j_idx+1):(j_idx+(D-1))]
    distance_mat[i,j] <- mat_dist(A, B, use_Riemann=use_Riemann)
    distance_mat[j,i] <- distance_mat[i,j]
  }
}

fit <- cmdscale(distance_mat, eig=TRUE, k=2) # k is the number of dim
cat("Lambda 1:",fit$eig[1],"\n")
cat("Lambda 2:",fit$eig[2],"\n")
cat("Lambda 3:",fit$eig[3],"\n")

plot_Sigma_ordination(fit, indiv_labels, "indiv", legend=FALSE)
plot_Sigma_ordination(fit, group_labels, "group")
