library(stray)
library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)
data(mallard_family) # this is a phyloseq object

# Just take vessel 1
ps <- prune_samples(sample_data(ps)$Vessel==1, ps)

# Just take hourly samples
ps <- prune_samples((sample_data(ps)$time > "2015-11-20 15:00:00 UTC") &
                      (sample_data(ps)$time < "2015-11-25 16:00:00 UTC"), ps)

# Order samples - to make plotting easy later
o <- order(sample_data(ps)$time) # get the order of a data.frame based on a column's value
otu_table(ps) <- otu_table(ps)[o,] # extract counts in this order; samples are rows!
sample_data(ps) <- sample_data(ps)[o,] # extract sample data in this order; samples are rows again -- phyloseq default?

# Extract Data / dimensions from Phyloseq object
Y <- t(as(otu_table(ps), "matrix")) # D x N
rownames(Y) <- taxa_names(ps)
D <- ntaxa(ps)
N <- nrow(sample_data(ps))

# X in hours
X <- as.numeric(sample_data(ps)$time)
X <- t((X-min(X)) / 3600) # 1 covariate (hour-number) over all samples

# Specify Priors
Gamma <- function(X) SE(X, sigma=5, rho=10) # Create partial function 
Theta <- function(X) matrix(0, D-1, ncol(X))
upsilon <- D-1+3
Xi <- matrix(.4, D-1, D-1)
diag(Xi) <- 1

# predict not just missing days but also forecast into future
X_predict <- t(1:(max(X))) # time point range, fill in any missing
predicted <- predict(fit.clr, X_predict, jitter=1) # predicts SAMPLES from the posterior (default = 2000)

family_names <- as(tax_table(ps)[,"Family"], "vector")
Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>% 
  gather_array(mean, coord, sample) %>% 
  mutate(time = X[1,sample], 
         coord = paste0("CLR(", family_names[coord],")"))

predicted_tidy <- gather_array(predicted, val, coord, sample, iter) %>% 
  mutate(time = X_predict[1,sample]) %>% # just adds time point ID
  filter(!is.na(val)) %>% # kill NAs
  group_by(time, coord) %>% # there are D x N = 1150 combos of time and coordinate (log ratio)
  summarise_posterior(val, na.rm=TRUE) %>% # gets standard quantiles for each coord x timepoint combo
  ungroup() %>% # just undoes this group x coord binding
  mutate(coord = paste0("CLR(", family_names[coord],")")) # replace coord with readable family name

ggplot(predicted_tidy, aes(x = time, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
  geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9)+
  geom_line(color="blue") +
  geom_point(data = Y_clr_tidy, alpha=0.5) +
  facet_wrap(~coord, scales="free_y") +
  theme_minimal()+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=45))










