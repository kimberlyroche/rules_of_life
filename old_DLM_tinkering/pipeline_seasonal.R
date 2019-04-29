source("include.R")

# what do the actual relative abundances look like?
if(FALSE) {
  glom_data <- load_glommed_data(level="species", replicates=TRUE)
  filtered <- filter_data(glom_data, count_threshold=3, sample_threshold=0.2, verbose=TRUE)
  
  ACA_data <- subset_samples(filtered, sname == "ACA")
  ACA_metadata <- read_metadata(ACA_data)
  
  n <- nsamples(ACA_data)
  p <- ntaxa(ACA_data)
  cat(paste("Data is",n,"samples by",p,"taxa!"))
  
  clr_data <- apply_clr(ACA_data) # 123 samples x 110 log ratios
  
  # plot 5 LR at random
  tax_idx <- sample(p)[sample(p)[1:5]]
  
  # row (sample ID) and column (lr taxon) order appears to be maintained but this
  # should be explicitly checked
  df <- gather_array(clr_data[,tax_idx], "log_rel_abundance", "sample_id", "taxon")
  df$sample_id <- factor(df$sample_id, levels=unique(df$sample_id))
  df$taxon <- factor(df$taxon, levels=unique(df$taxon))
  
  p <- ggplot(df, aes(x=sample_id, y=log_rel_abundance, group=taxon)) +
    geom_line(aes(color=taxon)) +
    geom_point(aes(color=taxon)) +
    theme_minimal() + 
    theme(legend.position="none")
  p
  #ggsave(paste("plots/ACA_clr_sample_plot.png",sep=""), plot=p, scale=2, width=10, height=4, units="in")
  
  # values in log ratio scale are approx. (-5, 5); exponentiated approx. (0.00674, 148)
  # remember these are ratios
}

# simulate log ratios with a Fourier-form rotational matrix G
# 5-component community, single harmonic

omega <- 2*pi/10 # period of 10 units of t

variance_scale <- 0.1
W <- diag(2) * variance_scale
#Xi <- diag(5) * variance_scale
#upsilon <- 10
Sigma <- diag(5) * variance_scale
gamma <- 1

F <- matrix(c(1, 0), 1, 2)
G <- matrix(c(cos(omega), -sin(omega), sin(omega), cos(omega)), 2, 2)

theta.t <- matrix(rnorm(10, 0, 2.5), 2, 5) # estimated from data above

T <- 100

etas <- matrix(0, T, 5)
for(t in 1:T) {
  #Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  theta.t <- G%*%theta.t + rmatrixnorm(matrix(0, 2, 5), W, Sigma)
  etas[t,] <- F%*%theta.t + rmvnorm(1, rep(0,5), gamma*Sigma)
}

df <- gather_array(etas, "log_rel_abundance", "time", "taxon")
df$taxon <- factor(df$taxon, levels=unique(df$taxon))

p <- ggplot(df, aes(x=time, y=log_rel_abundance, group=taxon)) +
  geom_line(aes(color=taxon)) +
  geom_point(aes(color=taxon)) +
  theme_minimal() + 
  theme(legend.position="none")
p
# amplitude explodes; this (see Pollard) is a known issue with these models at long time scales
# this happens because ???

# [TO DO] can we reason about how changes to Sigma and W affect the outcome?

# simulate log ratios with an oscillation in the observation
# 5-component community, single harmonic

omega <- 2*pi/10 # period of 10 units of t

W <- diag(2) * 0.01
#Xi <- diag(5) * variance_scale
#upsilon <- 10
Sigma <- diag(5) * 0.2
gamma <- 1

F.t <- matrix(0, T, 2)
for(t in 1:T) {
  F.t[t,] <- c(cos(omega*t), sin(omega*t))
}
G <- diag(2)

theta.t <- matrix(rnorm(10, 0, 0.1), 2, 5) # ?

T <- 100

etas <- matrix(0, T, 5)
for(t in 1:T) {
  #Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  theta.t <- G%*%theta.t + rmatrixnorm(matrix(0, 2, 5), W, Sigma)
  etas[t,] <- F%*%theta.t + rmvnorm(1, rep(0,5), gamma*Sigma)
}

df <- gather_array(etas, "log_rel_abundance", "time", "taxon")
df$taxon <- factor(df$taxon, levels=unique(df$taxon))

p <- ggplot(df, aes(x=time, y=log_rel_abundance, group=taxon)) +
  geom_line(aes(color=taxon)) +
  geom_point(aes(color=taxon)) +
  theme_minimal() + 
  theme(legend.position="none")
p

# [TO DO] probably useful to visualize this as proportions





