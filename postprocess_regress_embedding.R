source("include/R/GP.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
  stop("Arguments: (level)", call.=FALSE)
}
level <- args[1]

# regress embedded distance on these labels
#   group (factor)
#   matgroup
#   counts
#   density
#   mom (scratch: too many values)
#   dad (ditto)
#   momrank
#   drought
#   largegroup
#   momdied
#   competingsib
#   earlyadversity
#   birthrate_all
#   birthrate_surviving

df <- readRDS(paste0("output/plots/basset/",level,"/Sigma_ordination_centroids.rds"))
glom_data <- load_glommed_data(level=level, replicates=TRUE)

x_factors <- list(c("group", 1),
          c("matgroup", 1), 
          c("momrank", 1), 
          c("drought", 1), 
          c("largegroup", 1), 
          c("momdied", 1), 
          c("competingsib", 1), 
          c("earlyadversity", 1), 
          c("birthrate_all", 1), 
          c("birthrate_surviving", 1),
          c("counts", 0),
          c("density", 0))

for(x_factor in x_factors) {
  cat("Effect of",x_factor,":",x_factor[1],"\n")
  labeled_df <- na.omit(get_other_labels(df, glom_data, unique(df$labels), annotation=x_factor[1]))
  for(embed_dim in 1:6) {
    y <- labeled_df[,embed_dim]
    if(x_factor[2] > 0) {
      x <- as.factor(labeled_df[,7])
    } else {
      x <- as.numeric(labeled_df[,7])
    }
    var_expl <- summary(lm(y ~ x))$adj.r.squared
    if(var_expl > 0.2) {
      cat("\tR squared in dimension",embed_dim,":",round(var_expl, 3),"\n")
    }
  }
}

