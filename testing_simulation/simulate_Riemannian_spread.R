# this file simulates adding various amounts of noise to a random matrix, calculating the Riemannian distance,
# and visualizing what this does to a MDS embedding

relative_path <- ".."

library(Rcpp)
library(matrixsampling)
library(gridExtra)
library(ape)
library(Rtsne)

sourceCpp(file.path(relative_path,"include/cpp/Riemann_dist.cpp"))

D <- 20
N <- 20

all_samples <- matrix(NA, D, D*((N+1)*3))

plot <- F

upsilon <- D + 10
baseline <- rinvwishart(1, upsilon, diag(D)*(upsilon - D - 1))[,,1] # baseline

append_plot <- function(plist, sample_mat, label) {
  df <- gather_array(sample_mat, "value", "taxon1", "taxon2")
  p <- ggplot(df, aes(taxon1, taxon2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low="darkblue", high="darkred", name="correlation") +
    theme(legend.position = "none") +
    ggtitle(label)
  plist[[length(plist) + 1]] <- p
  return(plist)
}

plist <- list()
plist <- append_plot(plist, baseline, "mean")

upsilon <- D + 2 + 50
samples_close <- rinvwishart(N, upsilon, baseline*(upsilon - D - 1))
distances <- c()
all_samples[,1:D] <- baseline
for(i in 1:N) {
  distances <- c(distances, Riemann_dist_pair(baseline, samples_close[,,i]))
  if(plot) {
    plist <- append_plot(plist, samples_close[,,i], "close")
  }
  offset <- D
  left <- offset + (D*(i-1)) + 1
  right <- left + D - 1
  all_samples[,left:right] <- samples_close[,,i]
}
cat("Mean dist:",round(mean(distances), 3),"\n")

plist <- append_plot(plist, baseline, "mean")
upsilon <- D + 2 + 19
samples_med <- rinvwishart(N, upsilon, baseline*(upsilon - D - 1))
distances <- c()
all_samples[,((D*(N+1))+1):(D*(N+2))] <- baseline
for(i in 1:N) {
  distances <- c(distances, Riemann_dist_pair(baseline, samples_med[,,i]))
  if(plot) {
    plist <- append_plot(plist, samples_med[,,i], "med")
  }
  offset <- D*(N+2)
  left <- offset + (D*(i-1)) + 1
  right <- left + D - 1
  all_samples[,left:right] <- samples_med[,,i]
}
cat("Mean dist:",round(mean(distances), 3),"\n")

plist <- append_plot(plist, baseline, "mean")
upsilon <- D + 2 + 9
samples_far <- rinvwishart(N, upsilon, baseline*(upsilon - D - 1))
distances <- c()
all_samples[,((D*(2*N+2))+1):(D*(2*N+3))] <- baseline
for(i in 1:N) {
  distances <- c(distances, Riemann_dist_pair(baseline, samples_far[,,i]))
  if(plot) {
    plist <- append_plot(plist, samples_far[,,i], "far")
  }
  offset <- D*((N+1)*2 + 1)
  left <- offset + (D*(i-1)) + 1
  right <- left + D - 1
  all_samples[,left:right] <- samples_far[,,i]
}
cat("Mean dist:",round(mean(distances), 3),"\n")

if(plot) {
  do.call("grid.arrange", c(plist, nrow=3, ncol=6))
}

# perform the embedding

dist_mat <- matrix(NA, (3*(N+1)), (3*(N+1)))
for(i in 1:(3*(N+1))) {
  for(j in 1:(3*(N+1))) {
    if(i <= j) {
      dist_mat[i,j] <- Riemann_dist_pair(all_samples[,(D*(i-1)+1):(D*i)],all_samples[,(D*(j-1)+1):(D*j)])
    } else {
      dist_mat[i,j] <- dist_mat[j,i]
    }
  }
}

dist_mat_Euclid <- matrix(NA, (3*(N+1)), (3*(N+1)))
for(i in 1:(3*(N+1))) {
  for(j in 1:(3*(N+1))) {
    if(i <= j) {
      dist_mat_Euclid[i,j] <- dist(
        rbind(
          c(all_samples[,(D*(i-1)+1):(D*i)]),
          c(all_samples[,(D*(j-1)+1):(D*j)])
        )
      )
    } else {
      dist_mat_Euclid[i,j] <- dist_mat_Euclid[j,i]
    }
  }
}

# MDS
embedding <- cmdscale(dist_mat, k=2)
df <- data.frame(x=embedding[,1], y=embedding[,2],
                 labels=c("baseline", rep("close", N),
                          "baseline", rep("med", N),
                          "baseline", rep("far", N)))
ggplot(df, aes(x=x, y=y, color=labels)) +
  geom_point(size=2)

# MDS - Euclidean
embedding <- cmdscale(dist_mat_Euclid, k=2)
df <- data.frame(x=embedding[,1], y=embedding[,2],
                 labels=c("baseline", rep("close", N),
                          "baseline", rep("med", N),
                          "baseline", rep("far", N)))
ggplot(df, aes(x=x, y=y, color=labels)) +
  geom_point(size=2)








