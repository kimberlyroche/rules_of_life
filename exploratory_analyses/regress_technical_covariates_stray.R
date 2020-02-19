# this file uses stray::pibble to regress the count data from replicates on a bunch of variables
#
# we're interested in seeing if batch variables (flow_cell_lane, plate, library_pool, extract_dna_conc_ng) have any
#   discernable effects on taxonomical relative abundances; tl;dr -- good news, no

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

level <- "genus"
data <- load_and_filter("genus")

replicates <- subset_samples(data, sample_status==2) # replicates only
md <- sample_data(data)
n <- phyloseq::nsamples(data)
p <- ntaxa(data)
cat("Number samples:",n,"\n")
cat("Number taxa:",p,"\n")

# testing only
sample_idx <- sample(1:n)[1:50]
subset_ids <- sample_data(data)$sample_id[sample_idx]
# otherwise
#subset_ids <- sample_data(data)$sample_id
subsetted <- subset_samples(data, sample_id %in% subset_ids)
md <- sample_data(subsetted)

Y <- t(otu_table(subsetted)@.Data) # rows are taxa, columns are samples
rownames(Y) <- paste("taxon",seq(1:length(rownames(Y))))
X <- t(model.matrix(~as.factor(md$sname) +
                      as.factor(md$season) +
                      as.factor(md$flow_cell_lane) +
                      as.factor(md$plate) +
                      as.factor(md$library_pool) +
                      as.factor(round(md$extract_dna_conc_ng,digits=-1))))
png(file.path(relative_path,plot_dir,"technical_covariates_design_matrix.png"))
image(X)
dev.off()

new_rownames <- c(rownames(X)[1])
for(i in 2:length(rownames(X))) {
  new_rownames[i] <- str_replace(rownames(X)[i], 'as\\.factor\\(md\\$sname\\)', "sname:")
  new_rownames[i] <- str_replace(new_rownames[i], 'as\\.factor\\(md\\$season\\)', "season:")
  new_rownames[i] <- str_replace(new_rownames[i], 'as\\.factor\\(md\\$flow_cell_lane\\)', "flow_cell_lane:")
  new_rownames[i] <- str_replace(new_rownames[i], 'as\\.factor\\(md\\$plate\\)', "plate:")
  new_rownames[i] <- str_replace(new_rownames[i], 'as\\.factor\\(md\\$library_pool\\)', "library_pool:")
  new_rownames[i] <- str_replace(new_rownames[i], 'as\\.factor\\(round\\(md\\$extract_dna_conc_ng, digits = -1\\)\\)', "dna_conc_ng:")
}
rownames(X) <- new_rownames

D <- dim(Y)[1]
Q <- dim(X)[1]
N <- dim(Y)[2]
prior_obj <- default_ALR_prior(D)
upsilon <- prior_obj$upsilon
Theta <- matrix(0, D-1, Q)
Gamma <- diag(Q)
Xi <- prior_obj$Xi
cat("D x N x Q:",D,"x",N,"x",Q,"\n")
fit <- pibble(Y=Y, X=X, upsilon=upsilon, Theta=Theta, Gamma=Gamma, Xi=Xi, n_samples=2000)

#cat("Plotting posterior predictive fit...\n")
#ppc(fit) + ggplot2::coord_cartesian(ylim=c(0, 30000))
#ggsave(file.path(relative_path,plot_dir,"technical_covariates_pp.png"), scale=1)

fit_summary <- summary(fit, pars="Lambda")$Lambda
focus <- fit_summary[sign(fit_summary$p2.5) == sign(fit_summary$p97.5),]
focus <- unique(focus$coord)

# plot covariates in chunks
cat("Covariates are:\n")
print(rownames(X))

# these are sample subset specific!
plot(fit, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[c(2:5,45)])
ggsave(file.path(relative_path,plot_dir,"technical_covariates_individual_season.png",sep=""), scale=1.5)



