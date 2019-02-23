library(phyloseq)
library(ggplot2)

source("include.R")

translate_month <- function(x) {
  if(x == 2) { return("Feb") }
  #if(x == 3) { return("Mar") }
  if(x == 4) { return("Apr") }
  if(x == 5) { return("May") }
  if(x == 6) { return("Jun") }
  if(x == 7) { return("Jul") }
}

data <- readRDS("original_data/emp_baboon_pool_T_w_techReps.RDS")
md <- sample_data(data)
replicate_months <- unlist(lapply(md[md$sample_status==2,"collection_date"]@.Data, function(x) substr(x, 6, 7)))
reps_int <- as.numeric(replicate_months)

month_data <- data.frame(month=as.factor(unlist(lapply(reps_int, translate_month))))
levels(month_data$month) <- c("Feb", "Apr", "May", "Jun", "Jul")
p <- ggplot(month_data, aes(x=month)) + geom_bar()
ggsave("plots/replicate_frequency.png", plot=p, scale=1.5, height=4, width=4, units="in")

# see how much batch variables vary across these time points
batch_vars <- c("plate", "extract_dna_conc_ng", "flow_cell_lane", "library_pool")

replicates <- md[md$sample_status==2,]
replicates <- replicates[order(replicates$collection_date),]
replicates <- replicates[order(replicates$sname),]
write.table(t(replicates[,c("sample_id","collection_date","sname","sex","age","season","plate","extract_dna_conc_ng","library_pool","flow_cell_lane")]),
            file="temp.txt")

# reuse
replicates <- subset_samples(data, sample_status==2)
#glommed_reps_genus <- glom_counts(replicates, level="genus", NArm=FALSE)
glommed_reps_family <- glom_counts(replicates, level="family", NArm=FALSE)

individual <- "OMO"
indiv_samples <- subset_samples(glommed_reps_genus, sname==individual)
#indiv_samples <- subset_samples(glommed_reps_family, sname==individual)

# collapse all lowly expressed guys into one lump
taxa_sums <- colSums(otu_table(indiv_samples)@.Data)
taxa_IDs <- names(taxa_sums)
taxa_sums <- as.vector(taxa_sums)
merge_list <- taxa_IDs[which(taxa_sums < 25)]
merged_samples <- merge_taxa(indiv_samples, merge_list)
# append a rep # to each collection date, just for plotting
fake_dates <- sample_data(merged_samples)$collection_date
for(i in 1:length(fake_dates)) {
  fake_dates[i] <- paste(fake_dates[i],i,sep="_")
}
sample_data(merged_samples)$collection_date <- fake_dates
plot_timecourse(merged_samples, paste(individual,"_replicate_timecourse",sep=""), legend=T, legend_level="genus")
#plot_timecourse(merged_samples, paste(individual,"_replicate_timecourse",sep=""), legend=T, legend_level="family")
