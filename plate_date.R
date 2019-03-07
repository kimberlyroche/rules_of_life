source("include.R")

date_to_num <- function(date) {
  as.numeric(difftime(date, as.Date("2000-04-21"), units="days"))
}

load("glom_data_species_reps.RData")

md <- read_metadata(glom_data)
date <- sapply(md$collection_date, date_to_num)
plate <- md$plate
sname <- md$sname
data <- data.frame(date=date, plate=plate, sname=sname)
cat("R-squared date on plate:",summary(lm(plate~date, data))$r.squared,"\n")
cat("R-squared sname on plate:",summary(lm(plate~as.factor(sname), data))$r.squared,"\n")
p <- ggplot(data, aes(x=date, y=plate)) +
  geom_point() +
  theme_minimal()
p
ggsave("plots/plate_by_date_all.png", scale=1.5, width=4, height=4, units="in")

# what about correlation in replicates only?
reps <- subset_samples(glom_data, sample_status==2)
md <- read_metadata(reps)
date <- sapply(md$collection_date, date_to_num)
plate <- md$plate
sname <- md$sname
data <- data.frame(date=date, plate=plate, sname=sname)
cat("R-squared date on plate:",summary(lm(plate~date, data))$r.squared,"\n")
cat("R-squared sname on plate:",summary(lm(plate~sname, data))$r.squared,"\n")
p <- ggplot(data, aes(x=date, y=plate)) +
  geom_point() +
  theme_minimal()
p
ggsave("plots/plate_by_date_replicates.png", scale=1.5, width=4, height=4, units="in")
