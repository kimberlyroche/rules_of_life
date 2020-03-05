# this file simulates adding noise to a PD matrix and track changes in Riemannian distance to better understand
# "small" vs. "large" distances in this space
#
# here, we'll simulate linearly adding noise to a base matrix (A) by starting from random matrices A & B
# A^TA and (A+B)^T(A+B) will give base and noise-added PSD matrices respectively (important to watch that scale
# doesn't inadvertantly change here!)

relative_path <- ".."

library(ggplot2)
library(driver)

sourceCpp(file.path(relative_path,"include/cpp/Riemann_dist.cpp"))

# A is a baseline random matrix; we'll get a PSD guy via A^TA
# p is the proportion of noise to add
# scalar is just a total noise scale
add_noise_calc_dist <- function(A, p, scalar) {
  ATA <- t(A)%*%A # "baseline"
  m <- nrow(ATA)
  B <- matrix(rnorm(m*m, 0, scalar*sqrt(1/m)), m, m)
  AB <- (1-p)*A + p*B # titrate in noise
  ABTAB <- t(AB)%*%(AB)
  return(list(Riemannian=Riemann_dist_pair(ATA, ABTAB),
              Euclidean=dist(rbind(c(ATA), c(ABTAB))),
              samples=ABTAB))
}

# A is a baseline random matrix
# props are the proportions of noise to test (e.g. 0.1, 0.2, 0.3)
# no_samples is the number of sampled random noise additions to test
sweep_noise_proportions <- function(A, props, no_samples, scalar) {
  samples <- list()
  df <- data.frame(distance=c(), noise_prop=c(), measure=c()) # for plotting
  for(pp in 1:length(props)) {
    p <- props[pp]
    cat("Evaluating proportion:",round(p,3),"\n")
    Riemannian_distances <- c()
    Euclidean_distances <- c()
    samples[[pp]] <- array(NA, dim=c(m, m, no_samples))
    for(i in 1:no_samples) {
      noise_output <- add_noise_calc_dist(A, p, scalar)
      Riemannian_distances <- c(Riemannian_distances, noise_output$Riemannian)
      Euclidean_distances <- c(Euclidean_distances, noise_output$Euclidean)
      samples[[pp]][,,i] <- noise_output$sample
    }
    # grab 80% interval for the observed distances
    df <- rbind(df, data.frame(mean_distance=mean(Riemannian_distances),
                               min_distance=sort(Riemannian_distances)[round(no_samples*0.1)],
                               max_distance=sort(Riemannian_distances)[round(no_samples*0.9)],
                               noise_prop=p,
                               measure="Riemannian"))
    df <- rbind(df, data.frame(mean_distance=mean(Euclidean_distances),
                               min_distance=sort(Euclidean_distances)[round(no_samples*0.1)],
                               max_distance=sort(Euclidean_distances)[round(no_samples*0.9)],
                               noise_prop=p,
                               measure="Euclidean"))
  }
  return(list(samples=samples, df=df))
}

# --------------------------------------------------------------------------------------------------
# for demonstration purposes, let's spike in noise (to small matrices) and visualize this
# --------------------------------------------------------------------------------------------------

m <- 10
scalar <- 2
A <- matrix(rnorm(m*m, 0, scalar*sqrt(1/m)), m, m)
props <- c(0.05, 0.33, 0.99)
no_samples <- 20
sweep_output <- sweep_noise_proportions(A, props, no_samples, scalar)

# plot 5 samples at increasing noise levels
plist <- list()
df_ATA <- gather_array(cov2cor(t(A)%*%A), "value", "taxon1", "taxon2")
for(pp in 1:length(props)) {
  p <- ggplot(df_ATA, aes(taxon1, taxon2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low="darkblue", high="darkred", name="covariance") +
    ggtitle("baseline") +
    theme(legend.position = "none")
  plist[[length(plist)+1]] <- p
  for(i in 1:5) {
    df <- gather_array(cov2cor(sweep_output$samples[[pp]][,,i]), "value", "taxon1", "taxon2")
    p <- ggplot(df, aes(taxon1, taxon2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient2(low="darkblue", high="darkred", name="covariance") +
      ggtitle(paste0(round(props[pp]*100),"% noise")) +
      theme(legend.position = "none")
    plist[[length(plist)+1]] <- p
  }
}

p <- do.call("grid.arrange", c(plist, ncol=6))

# --------------------------------------------------------------------------------------------------
# test behavior of distance as f(noise) on a larger scale!
# --------------------------------------------------------------------------------------------------

m <- 50
n <- 20
A <- matrix(rnorm(m*m, 0, scalar*sqrt(1/m)), m, m)
props <- seq(0, 0.99, length.out=n)
no_samples <- 100
sweep_output <- sweep_noise_proportions(A, props, no_samples, scalar)

ggplot(sweep_output$df, aes(color=measure)) +
  geom_ribbon(aes(x=noise_prop, ymin=min_distance, ymax=max_distance), alpha=0.33) +
  geom_path(aes(x=noise_prop, y=mean_distance), size=1)

# --------------------------------------------------------------------------------------------------
# ordinate the output of the sweep so we can see how these distances behave;
# I don't know how instructive these really all
# --------------------------------------------------------------------------------------------------

include_sweeps <- c(1)
include_sweeps <- c(1, 2, 3)

all_samples <- t(A)%*%A
labels <- c("truth")
for(s in include_sweeps) {
  for(i in 1:n) {
    all_samples <- cbind(all_samples, sweep_output$samples[[s]][,,i])
    if(s == 1) {
      labels <- c(labels, "0.05 noise")
    } else if(s == 2) {
      labels <- c(labels, "0.33 noise")
    } else {
      labels <- c(labels, "0.99 noise")
    }
  }
}
d <- Riemann_dist_samples(all_samples, (length(include_sweeps)*n)+1, 1)

embedding <- cmdscale(d, k=2)
df <- data.frame(x1=embedding[,1], x2=embedding[,2],
                 type=labels)
fitted_point_sz <- 2
truth_point_size <- 4
p1 <- ggplot() +
  geom_point(data=df[df$type=="truth",], aes(x=x1, y=x2), size=(truth_point_size*1.3), shape=17) +
  geom_point(data=df[df$type!="truth",], aes(x=x1, y=x2, color=type), size=fitted_point_sz) +
  xlab("PCoA 1") +
  ylab("PCoA 2")
p1

# old: calculate distances over all proportions x samples; distance matrix with have this order
# note, this code is crazy inefficient
#
# -------------------------------------------------------------------------------------------------------
# | 0.1 samp 1 | ... | 0.1 samp N | 0.2 samp 1 | ... | 0.2 samp N | ... | 0.9 samp 1 | ... | 0.9 samp N |
# -------------------------------------------------------------------------------------------------------
# | 0.1 samp 2 | ...
# --------------

# calculate distances
# dist_mat_Riemannian <- matrix(NA, length(props)*no_samples, length(props)*no_samples)
# dist_mat_Euclidean <- matrix(NA, length(props)*no_samples, length(props)*no_samples)
# for(pp1 in 1:length(props)) {
#   p1 <- props[pp1]
#   offset1 <- no_samples*(pp1-1)
#   for(pp2 in 1:length(props)) {
#     p2 <- props[pp2]
#     cat(pp1,"x",pp2,"\n")
#     offset2 <- no_samples*(pp2-1)
#     for(i in 1:no_samples) {
#       for(j in 1:no_samples) {
#         dist_mat_Riemannian[(offset1+i),(offset2+j)] <- Riemann_dist_pair(samples[[pp1]][,,i], samples[[pp2]][,,j])
#         dist_mat_Euclidean[(offset1+i),(offset2+j)] <- dist(rbind(c(samples[[pp1]][,,i]), c(samples[[pp2]][,,j])))
#       }
#     }
#   }
# }

# match labels with structure diagrammed above
# labels <- c()
# for(pp in 1:length(props)) {
#   labels <- c(labels, rep(pp, no_samples))
# }

# embed Riemannian distances (horseshoe!)
# embedding <- cmdscale(dist_mat_Riemannian, k=2)
# df <- data.frame(x=embedding[,1], y=embedding[,2], label=as.factor(labels))
# ggplot(df) + 
#   geom_point(aes(x=x, y=y, color=label))

# embed Euclidean distances (another horseshoe!)
# embedding <- cmdscale(dist_mat_Euclidean, k=2)
# df <- data.frame(x=embedding[,1], y=embedding[,2], label=as.factor(labels))
# ggplot(df) +
#   geom_point(aes(x=x, y=y, color=label))

