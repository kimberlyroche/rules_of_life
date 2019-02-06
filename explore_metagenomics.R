library(driver)
library(dplyr)
library(ggplot2)

PC <- 0.65

# ============================================================================
# read 16S data
# ============================================================================

data_16S <- read.csv("original_data/Abundances_16S_unrarefied_OTU_table_AGG.csv",
                     header=TRUE,
                     stringsAsFactors=FALSE)

taxonomy <- data_16S[,"taxonomy"]
counts <- data_16S[, colnames(data_16S) != "taxonomy"]
rm(data_16S)
counts <- apply(counts, c(1,2), as.numeric)

hist_data <- gather_array(counts)
p <- hist_data %>% mutate(log_count=log(var+PC)) %>% ggplot() +
  geom_histogram(aes(x=log_count), binwidth=0.5)
p

# ============================================================================
# read shotgun metagenomics data for enzyme concentration
# ============================================================================

data_enzyme <- apply(read.csv("original_data/Abundances_KEGG_enzyme_orthologs.csv", 
                              header=TRUE, stringsAsFactors=FALSE), c(1,2), as.numeric)
data_enzyme <- data_enzyme[,colnames(data_enzyme) != "DIB"]

hist(as.vector(log(rowSums(data_enzyme)))) # what's the spread of the log rowsums?

levels <- as.vector(as.matrix(data_enzyme))
tiny_PC <- min(levels[levels != 0])*0.1
cat("Enzymes with total concentration below the pseudocount:",
    (sum(as.vector(rowSums(data_enzyme)) < tiny_PC)),"of",dim(data_enzyme)[1],"\n")

filter_idx <- apply(data_enzyme, 1, function(x) sum(x) > tiny_PC)
data_enzyme <- data_enzyme[as.vector(filter_idx),]

hist(as.vector(log(rowSums(data_enzyme))))

hist_data <- gather_array(data_enzyme)
p <- hist_data %>% mutate(log_conc=log(var+tiny_PC)) %>% ggplot() +
  geom_histogram(aes(x=log_conc), binwidth=0.5)
p

# ============================================================================
# read shotgun metagenomics data for KEGG pathways
# ============================================================================

data_pathways <- read.csv("original_data/Abundances_KEGG_pathways.csv", 
                          header=TRUE, stringsAsFactors=FALSE)
data_pathways <- data_pathways[,colnames(data_pathways) != "DIB"]

hist(as.vector(log(rowSums(data_pathways)))) # what's the spread of the low rowsums?

levels <- as.vector(as.matrix(data_pathways))
tiny_PC <- min(levels[levels != 0])*0.1
cat("Pathways with total concentration below the pseudocount:",
    (sum(as.vector(rowSums(data_pathways)) < tiny_PC)),"of",dim(data_pathways)[1],"\n")

filter_idx <- apply(data_pathways, 1, function(x) sum(x) > tiny_PC)
data_pathways <- data_pathways[as.vector(filter_idx),]

hist(as.vector(log(rowSums(data_pathways))))

hist_data <- gather_array(data_pathways)
p <- hist_data %>% mutate(log_conc=log(var+tiny_PC)) %>% ggplot() +
  geom_histogram(aes(x=log_conc), binwidth=0.5)
p

# ============================================================================
# read shotgun metagenomics data for KEGG modules
# ============================================================================

data_modules <- read.csv("original_data/Abundances_KEGG_modules.csv", 
                         header=TRUE, stringsAsFactors=FALSE)
data_modules <- data_modules[,colnames(data_modules) != "DIB"]

hist(as.vector(log(rowSums(data_modules)))) # what's the spread of the log rowsums?

levels <- as.vector(as.matrix(data_modules))
tiny_PC <- min(levels[levels != 0])*0.1
cat("Modules with total concentration below the pseudocount:",
    (sum(as.vector(rowSums(data_modules)) < tiny_PC)),"of",dim(data_modules)[1],"\n")

filter_idx <- apply(data_modules, 1, function(x) sum(x) > tiny_PC)
data_modules <- data_modules[as.vector(filter_idx),]

hist(as.vector(log(rowSums(data_modules))))

hist_data <- gather_array(data_modules)
p <- hist_data %>% mutate(log_conc=log(var+tinier_PC)) %>% ggplot() +
  geom_histogram(aes(x=log_conc), binwidth=0.5)
p

