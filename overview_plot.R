library(phyloseq)

date_to_num <- function(date) {
  as.numeric(difftime(date, as.Date("2000-04-21"), units="days"))
}

setwd("C:/Users/Kim/Documents/rules_of_life")
data <- readRDS("original_data/emp_baboon_pool_T_w_techReps.RDS")
metadata <- sample_data(data)

cat("Earliest sample:",min(metadata$collection_date),"\n")
cat("Latest sample:",max(metadata$collection_date),"\n")

#min_samples <- 50
min_samples <- 10

# get individuals with fewer than 50 samples and drop
individuals <- unique(metadata$sname)
samples <- list()
for(i in 1:length(individuals)) {
  indiv_samples <- subset_samples(metadata, (sname %in% c(individuals[i])))$collection_date
  indiv_samples <- as.vector(sapply(indiv_samples, date_to_num))
  if(length(indiv_samples) >= min_samples) {
    cat(individuals[i],"has",length(indiv_samples),"samples\n")
    samples[[individuals[i]]] <- indiv_samples
  }
}

plot.data <- data.frame(indiv=NULL, sample=NULL)
for(ind in names(samples)) {
  new.data <- data.frame(indiv=rep(ind, length(samples[[ind]])), sample=samples[[ind]])
  plot.data <- rbind(plot.data, new.data)
}

m <- ggplot(plot.data, aes(x=sample, y=indiv)) +
  geom_point() +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(margin=margin(t=10, r=0, b=0, l=10)))
#ggsave("C:/Users/kim/Desktop/ABRP_overview.png", scale=2, width=12, height=8, units="in")
ggsave("C:/Users/kim/Desktop/ABRP_overview.png", scale=2, width=4, height=8, units="in")
