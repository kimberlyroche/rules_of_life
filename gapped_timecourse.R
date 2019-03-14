library(RColorBrewer)

indiv <- "DUX"
indiv <- "DUI"
indiv <- "ACA"

data <- filter_data(data=glom_data, verbose=F)
data <- subset_samples(data, sname==indiv)
#data <- subset_samples(data, grepl("2002", collection_date))

cat("Min age:",min(sample_data(DUI)$age),"\n")
cat("Max age:",max(sample_data(DUI)$age),"\n")

save_filename <- paste(indiv,"_gapped_timecourse",sep="")
legend <- TRUE
legend_level <- "family"

p <- apply_proportion(data)
df <- psmelt(p)
df2 <- bind_cols(list(OTU=df$OTU, Sample=df$Sample, Abundance=df$Abundance))
  
# replace Sample ID's with their dates for readability
metadata <- sample_data(data)
for(i in 1:dim(df2)[1]) {
  df2$Sample[i] <- metadata[metadata$sample_id==df2$Sample[i],"collection_date"][[1]]
}

# make sure the dates are in order and fix the order by converting to factors
df3 <- df2[order(df2$Sample),]
gap.days <- 13
na.string <- ".N/A"
dates_present <- unique(df3$Sample)
for(d in 1:(length(dates_present)-1)) {
  diff <- as.Date(dates_present[d+1]) - as.Date(dates_present[d])
  next.date <- as.Date(dates_present[d])
  attr(diff, "units") <- "days"
  while(diff > gap.days) {
    next.date <- next.date + gap.days
    print(as.Date(dates_present[d]))
    print(next.date)
    df3 <- rbind(df3, list(OTU=na.string, Sample=as.character(next.date), Abundance=1))
    diff <- as.Date(dates_present[d+1]) - next.date
    print(diff)
    cat("\n")
  }
}
df3 <- df3[order(df3$Sample),]
df3$Sample <- factor(df3$Sample, levels=unique(df3$Sample))

# replace ASV sequences with their (abbrev.) taxonomy for readability
for(i in 1:dim(df3)[1]) {
  # show labels as order/family/genus
  # species is NA for all
  if(df3$OTU[i] != na.string) {
    if(legend_level == "species") {
      df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4:6]),collapse="/")
    } else if(legend_level == "genus") {
      df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4:5]),collapse="/")
    } else {
      df3$OTU[i] <- paste(as.vector(tax_table(data)[df3$OTU[i],4]),collapse="/")
    }
  }
}


categories <- unique(df3$OTU)
coul = brewer.pal(4, "Spectral")
coul = colorRampPalette(coul)(25)
coul[1] <- "#DDDDDD"
#palette(brewer.pal(n=length(categories), name="Spectral"))

p <- ggplot(df3, aes(x=Sample, y=Abundance, fill=OTU)) + 
  geom_bar(position="fill", stat="identity") +
  #scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=coul) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.text=element_text(size=8))
if(legend) {
  p <- p + theme(legend.position="bottom")
} else {
  p <- p + theme(legend.position="none")
}
ggsave(paste(save_filename,".png",sep=""), plot=p, scale=2, width=12, height=4, units="in")

















