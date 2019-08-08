source("include.R")

level <- "family"
glom_data <- load_glommed_data(level=level, replicates=TRUE)
data <- filter_data(glom_data, count_threshold=5, sample_threshold=0.33, collapse_level=level, verbose=TRUE)
tt <- tax_table(data)
rownames(tt) <- NULL
for(i in 1:(nrow(tt)-1)) {
  cat(tt[i,5],"/",tt[nrow(tt),5],"\n")
}

