library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)
library(scales)
library(coda)

# NEED TO FIGURE OUT WHY THIS DOES NOT WORK INSIDE A FUNCTION

#load("glom_data_g.RData")
#f2 <- filter_counts(glom_data, 3, 0.9)

individuals <- unique(read_metadata(f2)$sname)

lag.max <- 26
lag.sums <- numeric(lag.max)
lag.measured <- numeric(lag.max)
for(indiv in individuals) {
	counts <- subset_samples(f2, sname==indiv)
	log_ratios <- apply_ilr(counts)
	log_ratios <- t(apply(log_ratios, 1, function(x) x - mean(x)))
	cat(indiv,"has",dim(log_ratios)[2],"samples\n")
	md <- read_metadata(counts)
	# get distances between adjacent timepoints in weeks
	d1 <- as.Date(md$collection_date[1])
	week_d <- numeric(length(md$collection_date)-1)
	for(d in 2:length(md$collection_date)) {
		d2 <- as.Date(md$collection_date[d])
		week_d[d-1] <- round((d2-d1)/7) + 1
		d1 <- d2
	}
	# calculate max lag to use
	for(lag in 1:lag.max) {
		idx <- which(week_d==lag)
		if(length(idx) >= 1) {
			tot <- 0
			for(t in idx) {
				y.t <- as.vector(log_ratios[,t])
				y.tt <- sqrt(y.t%*%y.t)
				y.h <- as.vector(log_ratios[,t+1])
				y.hh <- sqrt(y.h%*%y.h)
				tot <- tot + (y.t%*%y.h)/(y.tt*y.hh)
			}
			lag.measured[lag] <- lag.measured[lag] + length(idx)
			lag.sums[lag] <- lag.sums[lag] + tot
		}
	}
}
for(lag	in 1:lag.max) {
	cat("Lag",lag,"=",(lag.sums[lag]/lag.measured[lag]),"\n")
	lag.sums[lag] <- lag.sums[lag]/lag.measured[lag]
}

p <- ggplot(as.data.frame(cbind(x=seq(1,lag.max), y=lag.sums)), aes(x=x, y=y)) +
	geom_line(linetype = "dashed") +
	scale_x_continuous(breaks=seq(0,lag.max,1)) +
	scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
	xlab("lag (weeks)") +
	ylab("ACF") +
	theme_minimal() +
	theme(axis.line.x=element_line(), axis.line.y=element_line()) +
	geom_hline(yintercept = 0)
ggsave("autocorrelation_all.png", plot=p, scale=2, width=4, height=3, units="in")
