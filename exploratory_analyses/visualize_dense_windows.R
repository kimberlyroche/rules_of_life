# this file plots the densely sampled windows identified by `find_dense_windows.R`

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

span <- "18mo"

data <- read.csv(file.path(relative_path,output_dir,paste0("windowed_sample_density_",span,".txt")), header=TRUE, sep="\t")

df <- data %>% 
  mutate(span=paste0(lower," -- ",upper)) %>%
  gather(indivs, density, dense10:dense30)
df$indivs <- as.factor(df$indivs)
levels(df$indivs) <- c(10, 15, 20, 25, 30)

p <- ggplot(df, aes(x=span, y=density, group=indivs)) +
	geom_line(aes(color=indivs)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
ggsave(file.path(relative_path,plot_dir,paste0("windows_",span,".png")), plot=p, dpi=150, units="in", height=6, width=20)


