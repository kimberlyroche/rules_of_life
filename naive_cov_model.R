library(driver)
library(mvtnorm)
library(MCMCpack)
library(ggplot2)

# global vars
pseudocount <- 0.65

n_samples <- 20 # number of posterior samples
interval_cutoff <- 0.8 # percent positive or negative required to be considered significantly non-zero

set.seed(1)

# https://stackoverflow.com/questions/25835643/replace-missing-values-with-column-mean
NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

# ---------------------------------------------------------------------------------------------------------------------
# parse paired sample data - 16S
# ---------------------------------------------------------------------------------------------------------------------

data_16S <- read.csv("original_data/Abundances_16S_unrarefied_OTU_table_AGG.csv", header=TRUE, stringsAsFactors=FALSE)
taxonomy <- data_16S[,"taxonomy"]
count_table <- data_16S[, colnames(data_16S) != "taxonomy"]
rm(data_16S)
count_table <- apply(count_table, c(1,2), as.numeric)

n <- dim(count_table)[2]
p1 <- dim(count_table)[1]

# ---------------------------------------------------------------------------------------------------------------------
# parse paired sample data - enzymes; OMIT FOR NOW, THERE'S ~9K
# ---------------------------------------------------------------------------------------------------------------------

# data_enzyme <- apply(read.csv("original_data/Abundances_KEGG_enzyme_orthologs.csv", 
#                               header=TRUE, stringsAsFactors=FALSE), c(1,2), as.numeric)
# data_enzyme <- data_enzyme[,colnames(data_enzyme) != "DIB"]
# data_enzyme <- log(data_enzyme)
# data_enzyme <- t(scale(t(data_enzyme)))

# p2 <- dim(data_enzyme)[1]

# ---------------------------------------------------------------------------------------------------------------------
# parse paired sample data - KEGG pathways
# ---------------------------------------------------------------------------------------------------------------------

data_pathways <- read.csv("original_data/Abundances_KEGG_pathways.csv", 
                          header=TRUE, stringsAsFactors=FALSE)
data_pathways <- data_pathways[,colnames(data_pathways) != "DIB"]
KEGG_labels <- rownames(data_pathways)

# remove all zero rows
zero_indices <- apply(data_pathways, 1, function(x) sum(x==0)==n)
data_pathways <- data_pathways[!(zero_indices),]

# log-transform and scale
levels <- as.vector(as.matrix(data_pathways))
tiny_PC <- min(levels[levels != 0])*0.1
data_pathways <- data_pathways + tiny_PC
data_pathways <- log(data_pathways+tiny_PC)
data_pathways <- t(scale(t(data_pathways)))

p3 <- dim(data_pathways)[1]

# # ---------------------------------------------------------------------------------------------------------------------
# # parse paired sample data - KEGG modules
# # ---------------------------------------------------------------------------------------------------------------------

data_modules <- read.csv("original_data/Abundances_KEGG_modules.csv", 
                         header=TRUE, stringsAsFactors=FALSE)
data_modules <- data_modules[,colnames(data_modules) != "DIB"]
KEGG_labels <- c(KEGG_labels, rownames(data_modules))

# remove all zero rows
zero_indices <- apply(data_modules, 1, function(x) sum(x==0)==n)
data_modules <- data_modules[!(zero_indices),]

# log-transform and scale
levels <- as.vector(as.matrix(data_modules))
tiny_PC <- min(levels[levels != 0])*0.1
data_modules <- data_modules + tiny_PC
data_modules <- log(data_modules+tiny_PC)
data_modules <- t(scale(t(data_modules)))

p4 <- dim(data_modules)[1]

# ---------------------------------------------------------------------------------------------------------------------
# fit model - CLR transform and plot correlation between features
# ---------------------------------------------------------------------------------------------------------------------

clr_counts <- count_table
for(subject in 1:n) {
  clr_counts[,subject] <- t(clr(as.vector(clr_counts[,subject])+pseudocount))
}

combined_data <- rbind(clr_counts, data_pathways, data_modules)

p <- p1 + p3 + p4

if(FALSE) {
  cat("Building covariance matrix...\n")
  cov_mat <- cov(t(combined_data))

  cat("Building correlation matrix...\n")
  cor_mat <- cov2cor(cov_mat)

  # plot correlation
  df <- gather_array(cor_mat)
  pl <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
    geom_tile() +
    xlab("") +
    ylab("") +
    theme_minimal() +
    guides(fill=guide_legend(title="correlation")) +
    ggtitle("Raw data correlation")
  ggsave("output/raw_data_correlation.png", plot=pl, scale=1.5, width=4, height=3, units="in")
}

# ---------------------------------------------------------------------------------------------------------------------
# fit model - set up weak prior over metabolites, microbes
# ---------------------------------------------------------------------------------------------------------------------

mu.0 <- matrix(0, nrow=p, ncol=1)
Lambda.0 <- diag(p)
#Lambda.0[(p1+1):p,(p1+1):p] <- diag(p2)*8 # this is the empirical variance of log-transformed OTU counts

nu.0 <- p + 2
kappa.0 <- 1

Sigmas <- array(dim=c(p, p, n_samples))
mus <- matrix(nrow=p, ncol=n_samples)

kappa.n <- kappa.0 + n
nu.n <- nu.0 + n

# ---------------------------------------------------------------------------------------------------------------------
# fit model - resample counts multinomial-Dirichlet
# ---------------------------------------------------------------------------------------------------------------------

for(s in 1:n_samples) {
  resampled_count_proportions <- count_table
  for(subject in 1:n) {
    alpha.j <- unlist(resampled_count_proportions[,subject]+pseudocount)
    resampled_count_proportions[,subject] <- t(rdirichlet(1, alpha.j))
  }
  # CLR transform
  clr_counts <- resampled_count_proportions
  rm(resampled_count_proportions)
  for(subject in 1:n) {
    clr_counts[,subject] <- t(clr(c(clr_counts[,subject])))
  }
  # combined dataset
  combined_data <- rbind(clr_counts, data_pathways, data_modules)

  # pre-calculate a few things
  ybar <- rowMeans(combined_data)
  x <- as.matrix(ybar - mu.0)
  dim(x) <- c(p, 1)
  S2 <- x%*%t(x)
  S <- matrix(0, nrow=p, ncol=p)
  for(j in 1:n) {
    x <- as.matrix(combined_data[,j] - ybar)
    dim(x) <- c(p, 1)
    S <- S + x%*%t(x)
  }
  Lambda.n <- Lambda.0 + S + ((kappa.0*n)/kappa.n)*S2

  # sample posterior Sigma (covariance matrix)
  Sigmas[,,s] <- riwish(nu.n, Lambda.n)

  # sample posterior theta (mean)
  mu.n <- (kappa.n*mu.0 + n*ybar)/kappa.n
  Sigma <- Lambda.n/(kappa.n*(nu.n-p+1))
  mus[,s] <- rmvt(1, sigma=Sigma, df=(nu.n-p+1))
}

# ---------------------------------------------------------------------------------------------------------------------
# select "significant" interactions
# ---------------------------------------------------------------------------------------------------------------------

# plot (average) correlation in resampled data
cor_mat <- cov2cor(apply(Sigmas, c(1,2), mean))
df <- gather_array(cor_mat)
pl <- ggplot(df, mapping=aes(x=dim_1, y=dim_2, fill=var)) +
            geom_tile() +
            xlab("") +
            ylab("") +
            theme_minimal() +
            guides(fill=guide_legend(title="correlation")) +
            ggtitle("Average correlation over samples of Sigma")
ggsave("output/avg_posterior_correlation.png", plot=pl, scale=1.5, width=4, height=3, units="in")

# convert covariance to correlation
cor_mat_samples <- array(dim=dim(Sigmas))
for (k in 1:n_samples) {
  cor_mat_samples[,,k] <- cov2cor(Sigmas[,,k])
}
cross_cor_samples <- cor_mat_samples[1:p1,(p1+1):p,]
rm(cor_mat_samples)

threshold_count <- round(interval_cutoff*n_samples)
filter_list <- apply(cross_cor_samples, c(1,2), function(x) abs(sum(sign(x))) > threshold_count)
cat("Percent above threshold:",sum(filter_list)/length(c(filter_list)),"\n")
# filter_list collapsed column-major where each column is the covariance of single 
# microbe with all KEGG modules/pathways
cov_list_df <- expand.grid(taxonomy, KEGG_labels)
cov_list <- apply(cov_list_df[, c("Var1","Var2")], 1, paste, collapse="#")
rm(cov_list_df)
rank_list <- list()
for (i in 1:p1) { # rows
  for (j in 1:(p3+p4)) { # columns
    idx <- (j-1)*p1 + i
    if(filter_list[i,j] == TRUE) {
      rank_list[cov_list[idx]] = mean(cross_cor_samples[i,j,])
    }
  }
}
rank_list <- rank_list[order(unlist(rank_list), decreasing=TRUE)]
save_file <- "output/ranked_list.RData"
save(data, Sigmas, rank_list, file=save_file)

# ---------------------------------------------------------------------------------------------------------------------
# format and return (nicely) the tagged interactions
# ---------------------------------------------------------------------------------------------------------------------

cat("Significant interactions:",length(rank_list),"\n")

metagenomics <- rbind(data_pathways, data_modules)

k <- 10
for(i in 1:k) {
  j <- length(rank_list)-i+1
  cat("RANK",i,"\n")
  cat("  Positive correlator\n")
  tax_ID <- strsplit(names(rank_list[i]), "#")[[1]][1]
  metagen_ID <- strsplit(names(rank_list[i]), "#")[[1]][2]
  cat("   ",tax_ID,"\n")
  cat("    ",metagen_ID,"\n")
  tax_idx <- which(taxonomy==tax_ID)
  metagen_idx <- which(KEGG_labels==metagen_ID)
  microbe_data <- clr_counts[tax_idx,]
  metagen_data <- metagenomics[metagen_idx,]
  df <- cbind(metagen_data=as.data.frame(metagen_data), microbe_data=microbe_data)
  sectioned_name <- paste(strsplit(tax_ID, ": ")[[1]][5:6], collapse="/")
  pl <- ggplot(df) +
    geom_point(aes(x=microbe_data, y=metagen_data)) +
    theme_minimal() +
    theme(text = element_text(size=8)) +
    xlab(paste(sectioned_name," (log counts)",sep="")) +
    ylab(paste(metagen_ID," (log conc)",sep=""))
  ggsave(paste("output/scatterplot_pos_",i,".png",sep=""), plot=pl, height=4, width=5)

  cat("  Negative correlator\n")
  tax_ID <- strsplit(names(rank_list[j]), "#")[[1]][1]
  metagen_ID <- strsplit(names(rank_list[j]), "#")[[1]][2]
  cat("    ",tax_ID,"\n")
  cat("    ",metagen_ID,"\n")
  tax_idx <- which(taxonomy==tax_ID)
  metagen_idx <- which(KEGG_labels==metagen_ID)
  microbe_data <- clr_counts[tax_idx,]
  metagen_data <- metagenomics[metagen_idx,]
  df <- cbind(metagen_data=as.data.frame(metagen_data), microbe_data=microbe_data)
  sectioned_name <- paste(strsplit(tax_ID, ": ")[[1]][5:6], collapse="/")
  pl <- ggplot(df) +
    geom_point(aes(x=microbe_data, y=metagen_data)) +
    theme_minimal() +
    theme(text = element_text(size=8)) +
    xlab(paste(sectioned_name," (log counts)",sep="")) +
    ylab(paste(metagen_ID," (log conc)",sep=""))
  ggsave(paste("output/scatterplot_neg_",i,".png",sep=""), plot=pl, height=4, width=5)
}










