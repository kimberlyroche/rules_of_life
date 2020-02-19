# this file renders a heatmap of where rows are posteriors of mean microbe x microbe correlation
# and columns are individual hosts
#
# this was EDA; untested for a while

relative_path <- ".."

library(RColorBrewer) 
library(gplots)

source(file.path(relative_path,"include/R/general.R"))

testing <- TRUE

level <- "genus"

# read all models
fns <- list.files(path=file.path(relative_path,model_dir,level))
models <- list()
ref_model <- NULL
chosen_model_idx <- 1:length(fns)
if(testing) {
  # just use 10 random individuals
  chosen_model_idx <- sample(length(fns))[1:10]
}
snames <- substr(fns, 1, 3)
R <- c()
for(i in chosen_model_idx) {
  cat("Reading model:",fns[i],"\n")
  R <- c(R, snames[i])
  temp <- readRDS(file.path(relative_path,model_dir,level,fns[i]))
  if(is.null(ref_model)) {
    ref_model <- temp
  }
  temp.fit <- to_clr(temp$fit)
  temp.Sigma <- temp.fit$Sigma
  for(k in 1:dim(temp.Sigma)[3]) {
    # convert to correlation
    temp.Sigma[,,k] <- cov2cor(temp.Sigma[,,k])
  }
  # take element-wise mean
  models[[i]] <- apply(temp.Sigma, c(1,2), mean)
}
cat("Models read.\n")

# vectorize correlation between all pairs
no_microbes <- nrow(models[[chosen_model_idx[1]]])
M <- matrix(NA, ((no_microbes^2)/2 - no_microbes/2), length(chosen_model_idx))
for(i in 1:length(chosen_model_idx)) {
  temp <- models[[chosen_model_idx[i]]]
  M[,i] <- c(temp[upper.tri(temp, diag=FALSE)])
}
# get shitty labels
L <- matrix(NA, no_microbes, no_microbes)
for(i in 1:no_microbes) {
  for(j in 1:no_microbes) {
    L[i,j] <- paste0(i," x ",j)
  }
}
L <- L[upper.tri(L, diag=FALSE)]

M_use <- t(M)
L_use <- L
R_use <- R
if(testing) {
  M_use <- t(M[1:200,])
  L_use <- L[1:200]
}

# render heatmap, clustering over rows (interactions) and columns (individuals)
hmcol <- colorRampPalette(brewer.pal(9, "Spectral"))(100)

cat("Rendering heatmap...\n")
png(file.path(relative_path,plot_dir,paste0(level,"_shared_relationships.png"), height=1000, width=3000)
h <- heatmap.2(M_use,
  main = "Microbe-microbe correlation", # heat map title
  labCol = L_use,                       # row labels (numbered interactions)
  labRow = R_use,                           # row labels (individuals)
  density.info = "none",                # turns off density plot inside color legend
  trace = "none",                       # turns off trace lines inside the heat map
  margins = c(4,4),                     # widens margins around plot
  col = hmcol,                          # use on color palette defined earlier
#  dendrogram="both")
  dendrogram="column")
eval(h$call)
dev.off()

