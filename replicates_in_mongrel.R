library(driver)
library(ggplot2)
library(devtools)
library(phyloseq)
library(stringr)

source("include.R")

level <- "genus"
filter <- 0.9

load(paste("glom_data_",level,"_reps.RData",sep=""))
replicates <- subset_samples(glom_data, sample_status==2)

# filter data down to family or genus
filtered <- filter_data(replicates, sample_threshold=filter, verbose=T)
md <- sample_data(filtered)
n <- phyloseq::nsamples(filtered)
p <- ntaxa(filtered)
cat("Number samples:",n,"\n")
cat("Number taxa:",p,"\n")

# does collection date predict plate?
dates <- numeric(n)
for(i in 1:n) {
  # date to number
  dates[i] <- as.numeric(as.POSIXct(md$collection_date[i], format="%Y-%m-%d"))
}
data <- data.frame(date=dates, plate=as.numeric(md$plate))
fit <- lm(plate ~ date, data)

devtools::load_all("C:/Users/kim/Documents/mongrel")

# for testing
#subset_ids <- sample_data(filtered)$sample_id[50:70]

subset_ids <- sample_data(filtered)$sample_id
subsetted <- subset_samples(filtered, sample_id %in% subset_ids)
md <- sample_data(subsetted)

Y <- t(otu_table(subsetted)@.Data) # rows are taxa, columns are samples
rownames(Y) <- paste("taxon",seq(1:length(rownames(Y))))
X <- t(model.matrix(~as.factor(md$sname) +
                      as.factor(md$season) +
                      as.factor(md$flow_cell_lane) +
                      as.factor(md$plate) +
                      as.factor(md$library_pool) +
                      as.factor(round(md$extract_dna_conc_ng,digits=-1))))
png(paste("design_",level,"_filter_",filter,".png",sep=""))
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

upsilon <- dim(Y)[1] + 2
Theta <- matrix(0, dim(Y)[1]-1, dim(X)[1])
Gamma <- diag(dim(X)[1])
Xi <- diag(dim(Y)[1]-1)
cat("N:",dim(Y)[2],"\n")
cat("D:",dim(Y)[1],"\n")
cat("Q:",dim(X)[1],"\n")
fit <- mongrel(Y=Y, X=X, upsilon=upsilon, 
                      Theta=Theta, Gamma=Gamma, Xi=Xi,
                      init=random_mongrel_init(Y), n_samples=2000)
cat("Saving fitted mongrel object...\n")
save(fit, file=paste("mongrelfit_",level,"_filter_",filter,".RData", sep=""))
dev.off()

cat("Plotting posterior-predictive fit...\n")
ppc(fit) + ggplot2::coord_cartesian(ylim=c(0, 30000))
ggsave(paste("ppc_",level,"_filter_",filter,".png",sep=""), scale=1)
dev.off()

fit_summary <- summary(fit, pars="Lambda")$Lambda
focus <- fit_summary[sign(fit_summary$p2.5) == sign(fit_summary$p97.5),]
focus <- unique(focus$coord)
plot(fit, par="Lambda", focus.coord = focus, focus.cov = rownames(X))
ggsave(paste("lambda_",level,"_filter_",filter,".png",sep=""), scale=2)
dev.off()











