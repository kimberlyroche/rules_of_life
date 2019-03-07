source("include.R")

load("glom_data_genus_reps.RData")
filtered <- filter_data(data=glom_data, sample_threshold=0.2)
md <- read_metadata(filtered)

limit <- 2000
limit <- nsamples(filtered)
sample_idx <- sample(nsamples(filtered))[1:limit]

plate <- md$plate[sample_idx]

covar <- c("season", "date", "grp", "sname")

for(covariate in covar) {
  if(covariate == "date") {
    cov_data <- unlist(lapply(md$collection_date, function(x) as.integer(format(as.Date(x), "%Y%m%d"))))[sample_idx]
  } else if(covariate == "season") {
    cov_data <- md$season[sample_idx]
  } else if(covariate == "grp") {
    cov_data <- md$grp[sample_idx]
  } else if(covariate == "sname") {
    cov_data <- md$sname[sample_idx]
  }

  data <- data.frame(plate=scale(plate), covariate=cov_data)

  # fit linear model; ask about percent variance explained
  if(covariate == "date") {
    res <- lm(plate~cov_data, data)
  } else {
    res <- lm(plate~factor(cov_data), data)
  }
  cat("Variance explained by",covariate,":",(1-var(res$residuals)),"\n")
  cat("Check with model:",summary(res)$r.squared,"\n")
  #print(res)

  if(covariate == "date") {
    p <- data %>% ggplot(aes(x=plate, y=covariate)) +
      geom_point(size=1) +
      theme_minimal() +
      xlab("plate") +
      ylab(covariate)
    ggsave(paste("testcov_",covariate,".png",sep=""), plot=p, scale=1.5, width=3, height=3, units="in")
  }
}

# how well does individual predict date anyway?
date <- unlist(lapply(md$collection_date, function(x) as.integer(format(as.Date(x), "%Y%m%d"))))[sample_idx]
date <- scale(date)
indiv <- md$sname[sample_idx]
data <- data.frame(x=indiv, y=date)
res <- lm(y~factor(x), data)
cat("Variance in collection_date explained by sname :",(1-var(res$residuals)),"\n")
#cat("Check with model:",summary(res)$r.squared,"\n")
