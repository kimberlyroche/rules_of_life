# this file moves a sliding window across samples and finds the number of individuals
# with >= K samples in windows of 12, 18, or 24 months

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

data <- readRDS(file.path(relative_path,data_dir,"filtered_family_5_20.rds"))
metadata <- sample_data(data)
hosts <- unique(metadata$sname)

span <- "24mo"

fileout <- file.path(relative_path,output_dir,paste0("windowed_sample_density_",span,".txt"))

cat("lower\tupper\tmean\tsd\tdense10\tdense15\tdense20\tdense25\tdense30",
    file=fileout,
    append=FALSE,
    sep="\n")

lower_date <- as.Date(min(metadata$collection_date))
max_date <- as.Date(max(metadata$collection_date))

while(lower_date < max_date) {
  if(span == "12mo") {
    days <- 365
  }
  if(span == "18mo") {
    days <- round(365*1.5)
  }
  if(span == "24mo") {
    days <- round(365*2)
  }
  upper_date <- lower_date + as.difftime(days, unit="days")
  sample_dist <- c()
  for(h in hosts) {
    indiv_samples <- subset_samples(data, sname==h)
    result = tryCatch({
      indiv_samples <- subset_samples(indiv_samples, (collection_date >= lower_date & collection_date <= upper_date))
      n <- phyloseq::nsamples(indiv_samples)
    }, error = function(e) {
      n <- 0
    })
    sample_dist <- c(sample_dist, n)
  }
  mean_val <- mean(sample_dist)
  sd_val <- sd(sample_dist)
  dense10 <- sum(sample_dist >= 10)
  dense15 <- sum(sample_dist >= 15)
  dense20 <- sum(sample_dist >= 20)
  dense25 <- sum(sample_dist >= 25)
  dense30 <- sum(sample_dist >= 30)
  cat(paste0(lower_date,"\t",upper_date,"\t",mean_val,"\t",sd_val,"\t",
             dense10,"\t",dense15,"\t",dense20,"\t",dense25,"\t",dense30),
    file=fileout,
    append=TRUE,
    sep="\n")
  lower_date <- lower_date + as.difftime(14, unit="days")
}

cat("Plotting results...\n")

# dumb
data <- read.csv(file.path(relative_path,output_dir,paste0("windowed_sample_density_",span,".txt")), header=TRUE, sep="\t")

df <- data %>%
  mutate(span=paste0(lower," -- ",upper)) %>%
  gather(indivs, density, dense10:dense30)
df$indivs <- as.factor(df$indivs)
levels(df$indivs) <- c(10, 15, 20, 25, 30)

p <- ggplot(df, aes(x=span, y=density, group=indivs)) +
        geom_line(aes(color=indivs)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
ggsave(file.path(plot_dir,paste0("windows_",span,".png")), plot=p, dpi=150, units="in", height=6, width=20)


