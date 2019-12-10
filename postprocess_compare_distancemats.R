library(phyloseq)
library(driver)
library(vegan)
library(ggplot2)

source("include/R/GP.R")
source("include/R/metagenomics.R")

level <- "family"

indiv_obj <- get_fitted_modellist_details(level=level)
hosts <- indiv_obj$individuals

cat("Building baselines...\n")
data <- readRDS(paste0("data/filtered_",level,"_5_20.rds"))
metadata <- sample_data(data)
tax <- tax_table(data)@.Data
indiv_mat <- NULL
indiv_labels <- c()
for(h in hosts) {
  cat("\thost:",h,"\n")
  indiv_labels <- c(indiv_labels, h)
  h_samples <- subset_samples(data, sname==h)
  #cat("\t# samples:",phyloseq::nsamples(h_samples),"\n")
  h_counts <- otu_table(h_samples)@.Data
  clr_counts <- driver::clr(h_counts + 0.5)
  mean_clr <- apply(clr_counts, 2, mean)
  if(is.null(indiv_mat)) {
    indiv_mat <- matrix(mean_clr, 1, length(mean_clr))
  } else {
    indiv_mat <- rbind(indiv_mat, mean_clr)
  }
}
rownames(indiv_mat) <- hosts

labels <- as.vector(tax[,5])
for(i in 1:length(labels)) {
  if(is.na(labels[i])) {
    for(j in 4:1) {
      if(!is.na(tax[i,j])) {
        if(j == 1) {
          labels[i] <- paste0("CLR(kingdom ",tax[i,j],")")
        }
        if(j == 2) {
          labels[i] <- paste0("CLR(phylum ",tax[i,j],")")
        }
        if(j == 3) {
          labels[i] <- paste0("CLR(class ",tax[i,j],")")
        }
        if(j == 4) {
          labels[i] <- paste0("CLR(order ",tax[i,j],")")
        }
        break
      }
    }
  } else {
    labels[i] <- paste0("CLR(family ",labels[i],")")
  }
}

indiv_prop <- driver::clrInv(indiv_mat)
colnames(indiv_prop) <- labels

exclude_low_prop <- as.vector(apply(indiv_prop, 2, function(x) mean(x) < 0.01))
indiv_prop <- indiv_prop[,!exclude_low_prop]

retained_list <- which(exclude_low_prop == FALSE)

df <- indiv_prop %>% gather_array(value, sname, taxon)
df$sname <- as.factor(df$sname)
levels(df$sname) <- rownames(indiv_prop)
df$taxon <- as.factor(df$taxon)
levels(df$taxon) <- colnames(indiv_prop)

categories <- unique(df$taxon)
coul = brewer.pal(4, "Spectral")
coul = colorRampPalette(coul)(length(unique(df$taxon)))
p <- ggplot(df, aes(x=sname, y=value, fill=taxon)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=coul) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.position="bottom")
ggsave(paste0(plot_dir,"average_baselines.png"), plot=p, dpi=100, width=16, height=5, units="in")

# sanity check these average
check_hosts <- c("AYU", "LIZ", "ZAI")
for(ch in check_hosts) {
  indiv_data <- subset_samples(data, sname == ch)
  indiv_prop <- otu_table(indiv_data)@.Data
  for(i in 1:nrow(indiv_prop)) {
    indiv_prop[i,] <- indiv_prop[i,]/sum(indiv_prop[i,])
  }
  colnames(indiv_prop) <- labels
  indiv_prop <- indiv_prop[,!exclude_low_prop]

  df <- indiv_prop %>% gather_array(value, sname, taxon)
  df$sname <- as.factor(df$sname)
  levels(df$sname) <- rownames(indiv_prop)
  df$taxon <- as.factor(df$taxon)
  levels(df$taxon) <- colnames(indiv_prop)

  categories <- unique(df$taxon)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df$taxon)))
  p <- ggplot(df, aes(x=sname, y=value, fill=taxon)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.position="bottom")
  p <- p + theme(legend.position="bottom")
  ggsave(paste0(plot_dir,"average_baseline_",ch,".png"), plot=p, dpi=100, width=10, height=8, units="in")
}

if(!file.exists(paste0(output_dir,"indiv_dist_1.rds")) | !file.exists(paste0(output_dir,"indiv_dist_2.rds"))) {

  # --------------------------------------------------------------------------------------
  #   fit (Euclidean/Aitchison) individual distances
  # --------------------------------------------------------------------------------------

  indiv_dist.1 <- dist(indiv_mat)
  attr(indiv_dist.1, "Labels") <- hosts
  saveRDS(indiv_dist.1, paste0(output_dir,"indiv_dist_1.rds"))

  # --------------------------------------------------------------------------------------
  #   fit Riemannian distance over individuals from MAP estimates of dynamics
  # --------------------------------------------------------------------------------------

  cat("Calculating DYNAMICS distances...\n")
  indiv_dist.2 <- calc_posterior_distances(level, which_measure="Sigma", indiv=NULL)
  indiv_dist.2 <- as.dist(indiv_dist.2$distance_mat)
  attr(indiv_dist.2, "Labels") <- hosts
  saveRDS(indiv_dist.2, paste0(output_dir,"indiv_dist_2.rds"))

} else {

  indiv_dist.1 <- readRDS(paste0(output_dir,"indiv_dist_1.rds"))
  indiv_dist.2 <- readRDS(paste0(output_dir,"indiv_dist_2.rds"))

}

which <- 1

if(which == 1) {
  dimx <- dim(as.matrix(indiv_dist.1))
  null_points <- matrix(rnorm(dimx[1]*ntaxa(data)), dimx[1], ntaxa(data))
  indiv_dist.null <- dist(null_points)
  # pure noise
  mantel(indiv_dist.null, indiv_dist.2, method="pearson", permutations=999)
  t1 <- as.matrix(indiv_dist.null)
  t1 <- t1[upper.tri(t1, diag=F)]
  t2 <- as.matrix(indiv_dist.2)
  t2 <- t2[upper.tri(t2, diag=F)]
  df <- data.frame(x=t1, y=t2)
  p <- ggplot(df, aes(x=x, y=y)) +
    geom_point()
  ggsave(paste0(plot_dir,"mantel_noise.png"), plot=p, dpi=100, units="in", height=8, width=8)
  
  # apparent correlation
  mantel(indiv_dist.1, indiv_dist.2, method="pearson", permutations=999)
  t1 <- as.matrix(indiv_dist.1)
  t1 <- t1[upper.tri(t1, diag=F)]
  df <- data.frame(x=t1, y=t2)
  p <- ggplot(df, aes(x=x, y=y)) +
    geom_point()
  ggsave(paste0(plot_dir,"mantel_test.png"), plot=p, dpi=100, units="in", height=8, width=8)
}

if(which == 2) {
   h1 <- hclust(indiv_dist.1)
   h2 <- hclust(indiv_dist.2)

   order1 <- h1$labels[h1$order]
   order2 <- h2$labels[h2$order]

  df.baseline.1 <- data.frame(host=c(), pair=c(), value=c())
  df.baseline.2 <- data.frame(host=c(), pair=c(), value=c())
  df.dynamics.1 <- data.frame(host=c(), pair=c(), value=c())
  df.dynamics.2 <- data.frame(host=c(), pair=c(), value=c())
  for(h in 1:length(order1)) {
    h.1 <- order1[h]
    h.2 <- order2[h]
    # baseline, ordered by dist(baseline)
    temp_vec <- indiv_mat[h.1,]
    df.baseline.1 <- rbind(df.baseline.1, data.frame(host=h.1, coord=1:length(temp_vec), value=temp_vec))
    # baseline, ordered by dist(dynamics)
    temp_vec <- indiv_mat[h.2,]
    df.baseline.2 <- rbind(df.baseline.2, data.frame(host=h.2, coord=1:length(temp_vec), value=temp_vec))
    # dynamics, ordered by dist(baseline)
    temp <- readRDS(paste0("output/model_fits/family/",h.1,"_bassetfit.rds"))$fit$Sigma
    temp_vec <- c(temp[upper.tri(temp, diag=F)])
    df.dynamics.1 <- rbind(df.dynamics.1, data.frame(host=h.1, pair=1:length(temp_vec), value=temp_vec))
    # dynamics, ordered by dist(dynamics)
    temp <- readRDS(paste0("output/model_fits/family/",h.2,"_bassetfit.rds"))$fit$Sigma
    temp_vec <- c(temp[upper.tri(temp, diag=F)])
    df.dynamics.2 <- rbind(df.dynamics.2, data.frame(host=h.2, pair=1:length(temp_vec), value=temp_vec))
  }

  p1 <- ggplot(df.baseline.1, aes(coord, host)) +
      geom_tile(aes(fill=value), colour="white") +
      scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
  ggsave(paste0(plot_dir,"hclust_ordered_baseline_by_baseline.png"), plot=p1, dpi=100, units="in", height=12, width=20)

  p2 <- ggplot(df.baseline.2, aes(coord, host)) +
      geom_tile(aes(fill=value), colour="white") +
      scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
  ggsave(paste0(plot_dir,"hclust_ordered_baseline_by_dynamics.png"), plot=p2, dpi=100, units="in", height=12, width=20)

  p3 <- ggplot(df.dynamics.1, aes(pair, host)) +
      geom_tile(aes(fill=value), colour="white") +
      scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
  ggsave(paste0(plot_dir,"hclust_ordered_dynamics_by_baseline.png"), plot=p3, dpi=100, units="in", height=12, width=20)

  p4 <- ggplot(df.dynamics.2, aes(pair, host)) +
      geom_tile(aes(fill=value), colour="white") +
      scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
  ggsave(paste0(plot_dir,"hclust_ordered_dynamics_by_dynamics.png"), plot=p4, dpi=100, units="in", height=12, width=20)
}
