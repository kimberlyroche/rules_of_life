#source("include.R")
#filtered <- filter_data(sample_threshold=0.9)
#md <- read_metadata(filtered)

limit <- 5000
sample_idx <- sample(nsamples(filtered))[1:limit]

y <- as.factor(md$plate)[sample_idx]
y2 <- md$plate[sample_idx]

#x <- unlist(lapply(md$collection_date, function(x) as.integer(format(as.Date(x), "%Y%m%d"))))[sample_idx]
x <- as.factor(md$sname)[sample_idx]

png("test.png")
plot(x, y, xlab="date", ylab="plate")
dev.off()

#cat("Spearman cor:",cor(x, y2, method="spearman"),"\n")
