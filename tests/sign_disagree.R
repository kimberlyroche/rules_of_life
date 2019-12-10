library(stray)

b1 <- "LYE"
b2 <- "NIN"

fit1 <- readRDS(paste0("output/model_fits/genus_SAVED/",b1,"_bassetfit.rds"))
fit2 <- readRDS(paste0("output/model_fits/genus_SAVED/",b2,"_bassetfit.rds"))

fit1.clr <- to_clr(fit1$fit)
fit2.clr <- to_clr(fit2$fit)

S1 <- apply(fit1.clr$Sigma, c(1,2), mean)
S2 <- apply(fit2.clr$Sigma, c(1,2), mean)

S1.tri <- S1[upper.tri(S1, diag=TRUE)]
S2.tri <- S2[upper.tri(S2, diag=TRUE)]

S1.pos <- (S1 >= 0)
S2.pos <- (S2 >= 0)
percent_disagree <- (1 - sum(S1.pos == S2.pos)/length(S1.pos))
cat(paste0("Percent disagreement between ",b1," and ",b2,": ",round(percent_disagree*100),"\n"))


