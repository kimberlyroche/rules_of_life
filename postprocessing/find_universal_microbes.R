# this file represents one attempt to identify "universal" microbes (those will a set of associations to other
# taxa that is preserved [the set of relationships is preserved, that is] across hosts)
#
# note: need to re-run this and decide whether this approach is worth keeping

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

level <- "family"
# pull all fitted models
models <- get_fitted_modellist_details(level=level)

Sigmas <- list()

for(i in 1:length(models$individuals)) {
  host <- models$individuals[i]
  cat("Processing",host,"\n")
  fit <- readRDS(models$model_list[i])$fit
  dim(fit$Eta) <- c(dim(fit$Eta), 1)
  dim(fit$Lambda) <- c(dim(fit$Lambda), 1)
  dim(fit$Sigma) <- c(dim(fit$Sigma), 1)
  fit.clr <- to_clr(fit)
  Sigmas[[host]] <- fit.clr$Sigma[,,1]
}

ref <- readRDS(file.path(relative_path,data_dir,"filtered_family_5_20.rds"))
tax <- tax_table(ref)@.Data
labels <- as.vector(tax[,5])
for(i in 1:length(labels)) {
  if(i == 42) {
    # this is "Other" category for family-level data filtered at 5-counts in 20% or more samples
    labels[i] <- "CLR(Other)"
  } else {
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
}

# columns are taxa
ref.clr <- driver::clr(otu_table(ref)@.Data + 0.5)
mean.clr <- apply(ref.clr, 2, mean)

hosts <- models$individuals
correlations <- data.frame(taxon=c(), pair=c(), value=c(), meanclr=c())
for(i in 1:models$D) {
  for(h1 in 1:(length(hosts)-1)) {
    for(h2 in (h1+1):length(hosts)) {
      # unique pair of individuals
      host1 <- hosts[h1]
      host2 <- hosts[h2]
      cat("Comparing",host1,"and",host2,"on taxon",labels[i],"(",i,")...\n")
      correlations <- rbind(correlations,
                            data.frame(taxon=labels[i],
                                       pair=paste0(host1,"x",host2),
                                       value=cor(Sigmas[[host1]][i,], Sigmas[[host2]][i,]),
                                       meanclr=mean.clr[i]))
    }
  }
}

saveRDS(Sigmas, file.path(relative_path,output_dir,"Sigmas.rds"))
saveRDS(correlations, file.path(relative_path,output_dir,"correlations.rds"))

p <- ggplot(correlations, aes(value)) +
     geom_density(aes(color=meanclr)) +
     facet_wrap(~ taxon, ncol=10) +
     scale_color_gradient(low='blue', high='red')
ggsave(file.path(relative_path,plot_dir,"universal_microbes.png"), plot=p, dpi=100, units="in", height=10, width=20)

mean_ranked <- as.data.frame(correlations %>% group_by(taxon) %>% summarize(mean_value = mean(value)))
head(mean_ranked)

med_ranked <- as.data.frame(correlations %>% group_by(taxon) %>% summarize(median_value = median(value)))
head(med_ranked)

med_ranked[order(med_ranked$median_value),]


