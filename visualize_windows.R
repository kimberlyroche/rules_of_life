library(ggplot2)
library(tidyverse)

span <- "18mo"

data <- read.csv(paste0(output_dir,"windowed_sample_density_",span,".txt"), header=TRUE, sep="\t")

df <- data %>% 
  mutate(span=paste0(lower," -- ",upper)) %>%
  gather(indivs, density, dense10:dense30)
df$indivs <- as.factor(df$indivs)
levels(df$indivs) <- c(10, 15, 20, 25, 30)

p <- ggplot(df, aes(x=span, y=density, group=indivs)) +
	geom_line(aes(color=indivs)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
ggsave(paste0(plot_dir,"windows_",span,".png"), plot=p, dpi=150, units="in", height=6, width=20)


