library(ggplot2)

Sigmas <- readRDS(paste0(output_dir,"Sigmas.rds"))
hosts <- names(Sigmas)
D <- nrow(Sigmas[[hosts[1]]])

m1 <- 35
m1_label <- "family_Helicobacteraceae"
m2 <- 6
m2_label <- "order_Clostridiales"

df.m1 <- data.frame(coord=c(), host=c(), value=c())
df.m2 <- data.frame(coord=c(), host=c(), value=c())
for(h in hosts) {
	row <- Sigmas[[h]][m1,]
	df.m1 <- rbind(df.m1, data.frame(coord=1:D, host=rep(h, D), value=row))
	row <- Sigmas[[h]][m2,]
	df.m2 <- rbind(df.m2, data.frame(coord=1:D, host=rep(h, D), value=row))
}

p1 <- ggplot(df.m1, aes(coord, host)) +
  geom_tile(aes(fill=value), colour="white") +
  scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
ggsave(paste0(plot_dir,m1_label,".png"), plot=p1, dpi=100, units="in", height=14, width=20)

p2 <- ggplot(df.m2, aes(coord, host)) +
  geom_tile(aes(fill=value), colour="white") +
  scale_fill_gradient2(low="darkblue", mid="white", high="darkred")
ggsave(paste0(plot_dir,m2_label,".png"), plot=p2, dpi=100, units="in", height=14, width=20)

