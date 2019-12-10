library(stray)
library(phyloseq)
library(RColorBrewer) 
library(gplots)
library(isa2)

source("include/R/general.R")

# render a heatmap of:
#   rows: posterior mean microbe x microbe correlation
#   cols: individual

testing <- TRUE

level <- "genus_SAVED"

# read all models
fns <- list.files(path=paste0(model_dir, level))
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
  temp <- readRDS(paste0(model_dir, level, "/", fns[i]))
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

# bicluster via ISA
# via manual: https://cran.r-project.org/web/packages/isa2/isa2.pdf

## Just to get always the same result
set.seed(24)
row.seeds <- generate.seeds(length=nrow(M_use), count=round(nrow(M_use/5)))
col.seeds <- generate.seeds(length=ncol(M_use), count=round(ncol(M_use/5)))
normed.data <- isa.normalize(M_use)
isaresult <- isa.iterate(normed.data, thr.row=1, thr.col=1, row.seeds=row.seeds, col.seeds=col.seeds)

png("test.png")
image(isaresult$rows%*%t(isaresult$columns))
dev.off()

quit()

# render heatmap, clustering over rows (interactions) and columns (individuals)
hmcol <- colorRampPalette(brewer.pal(9, "Spectral"))(100)

if(FALSE) {
cat("Rendering heatmap...\n")
png("test.png", height=1000, width=3000)
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
}

# -----------------------------------------------------
# order individuals/rows by (1) group (2) sample number
# -----------------------------------------------------

if(FALSE) {
group <- c()
sample_number <- c()
data <- load_and_filter(level="phylum")
for(indiv in snames[chosen_model_idx]) {
  indiv <<- indiv
  indiv_subset <- subset_samples(data, sname==indiv)
  grp <- sample_data(indiv_subset)$grp
  ux <- unique(grp)
  group <- c(group, ux[which.max(tabulate(match(grp, ux)))][[1]])
  sample_number <- c(sample_number, phyloseq::nsamples(indiv_subset))
}

if(testing) {
  M_use <- t(M[1:200,])
  L_use <- L[1:200]
}
M_use <- M_use[order(sample_number),]
R_use <- paste0(R[order(sample_number)]," ",sample_number[order(sample_number)])

cat("Rendering sample-number-ordered heatmap...\n")
png("test2.png", height=1000, width=3000)
h <- heatmap.2(M_use,
  main = "Microbe-microbe correlation", # heat map title
  labCol = L_use,                       # row labels (numbered interactions)
  labRow = R_use,                           # row labels (individuals)
  density.info = "none",                # turns off density plot inside color legend
  trace = "none",                       # turns off trace lines inside the heat map
  margins = c(4,6),                     # widens margins around plot
  col = hmcol,                          # use on color palette defined earlier
#  dendrogram="both")
  Rowv=NULL,
  dendrogram="column")
eval(h$call)
dev.off()

if(testing) {
  M_use <- t(M[1:200,])
  L_use <- L[1:200]
}
M_use <- M_use[order(group),]
R_use <- paste0(R[order(group)]," ",group[order(group)])

cat("Rendering group-ordered heatmap...\n")
png("test3.png", height=1000, width=3000)
h <- heatmap.2(M_use,
  main = "Microbe-microbe correlation", # heat map title
  labCol = L_use,                       # row labels (numbered interactions)
  labRow = R_use,                           # row labels (individuals)
  density.info = "none",                # turns off density plot inside color legend
  trace = "none",                       # turns off trace lines inside the heat map
  margins = c(4,6),                     # widens margins around plot
  col = hmcol,                          # use on color palette defined earlier
#  dendrogram="both")
  Rowv=NULL,
  dendrogram="column")
eval(h$call)
dev.off()
}
