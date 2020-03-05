# this file represents one attempt to identify "universal" microbes (those will a set of associations to other
# taxa that is preserved [the set of relationships is preserved, that is] across hosts)
#
# note: need to re-run this and decide whether this approach is worth keeping

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

level <- "family" # TBD: repeat these results at the level of genera!
level_no <- 5

level <- "genus"
level_no <- 6

# pull all fitted models
models <- get_fitted_modellist_details(level=level, MAP=TRUE)

Sigmas.clr <- list()
Sigmas.alr <- list()

for(i in 1:length(models$individuals)) {
  host <- models$individuals[i]
  cat("Processing",host,"\n")
  fit <- readRDS(models$model_list[i])$fit
  dim(fit$Eta) <- c(dim(fit$Eta), 1)
  dim(fit$Lambda) <- c(dim(fit$Lambda), 1)
  dim(fit$Sigma) <- c(dim(fit$Sigma), 1)
  Sigmas.alr[[host]] <- fit$Sigma[,,1]
  fit.clr <- to_clr(fit)
  Sigmas.clr[[host]] <- fit.clr$Sigma[,,1]
}

# ======================================================================================================
#   LOOK FOR STRONGEST INTERACTIONS
#   we'll sanity check by plotting predictions in the ALR, then report in the CLR
# ======================================================================================================

# hard-coding to family level here
# if string `alr_ref` is provided, label.alr will be returned as
#   ALR(x / y)
# otherwise, label.alr will be returned as
#   x
# label.clr will always be returned as
#   CLR(x)
walk_up_NA_label <- function(i, tax, level_no, alr_ref=NULL) {
  label.clr <- NA
  label.alr <- NA
  for(j in level_no:1) {
    if(!is.na(tax[i,j])) {
      if(j == 1) {
        label.clr <- paste0("CLR(kingdom ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(kingdom ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("kingdom ",tax[i,j])
        }
      }
      if(j == 2) {
        label.clr <- paste0("CLR(phylum ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(phylum ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("phylum ",tax[i,j])
        }
      }
      if(j == 3) {
        label.clr <- paste0("CLR(class ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(class ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("class ",tax[i,j])
        }
      }
      if(j == 4) {
        label.clr <- paste0("CLR(order ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(order ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("order ",tax[i,j])
        }
      }
      if(j == 5) {
        label.clr <- paste0("CLR(family ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(family ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("family ",tax[i,j])
        }
      }
      if(j == 6) {
        label.clr <- paste0("CLR(genus ",tax[i,j],")")
        if(!is.null(alr_ref)) {
          label.alr <- paste0("ALR(genus ",tax[i,j]," / ",alr_ref,")")
        } else {
          label.alr <- paste0("genus ",tax[i,j])
        }
      }
      break
    }
  }
  return(list(clr=label.clr, alr=label.alr))
}

# build some readable labels
ref <- readRDS(file.path(relative_path,data_dir,paste0("filtered_",level,"_5_20.rds")))
tax <- tax_table(ref)@.Data
labels.clr <- as.vector(tax[,level_no]) # family
labels.alr <- as.vector(tax[,level_no])
alr.ref <- walk_up_NA_label(nrow(tax), tax, level_no)$alr
for(i in 1:length(labels.clr)) {
  if(i == 42) {
    # this is "Other" category for family-level data filtered at 5-counts in 20% or more samples
    labels.clr[i] <- "CLR(Other)"
    labels.alr[i] <- paste0("ALR(Other / ",alr.ref,")")
  } else {
    pieces <- walk_up_NA_label(i, tax, level_no, alr_ref=alr.ref)
    labels.clr[i] <- pieces$clr
    labels.alr[i] <- pieces$alr
  }
}
labels.clr

# ======================================================================================================
#   collect pairwise correlations and labels associated with them
# ======================================================================================================

# build a map from which we can translate an interaction number with the original pair of taxa
label_pairs.clr <- matrix(NA, models$D, models$D)
label_pairs.alr <- matrix(NA, models$D-1, models$D-1)
for(i in 1:models$D) {
  for(j in 1:models$D) {
    if(i < j) {
      label_pairs.clr[i,j] <- paste0(i,"_",j)
      if(i < models$D & j < models$D) {
        label_pairs.alr[i,j] <- paste0(i,"_",j)
      }
    }
  }
}

# plot all in heatmap as hosts x pairs of microbial interaction
# `.tri` here indicates we're just keeping the upper triangular of the interaction matrix, i.e.
#   the unique interactions
interaction_pairs.clr <- matrix(NA, length(models$individuals), models$D^2)
interaction_pairs.alr <- matrix(NA, length(models$individuals), (models$D-1)^2)
interaction_pairs.clr.tri <- matrix(NA, length(models$individuals), (models$D^2)/2 - models$D/2)
interaction_pairs.alr.tri <- matrix(NA, length(models$individuals), ((models$D-1)^2)/2 - (models$D-1)/2)
for(i in 1:length(models$individuals)) {
  host <- models$individuals[i]
  # convert this host's MAP covariance to correlation
  corr_mat.clr <- cov2cor(Sigmas.clr[[host]])
  corr_mat.alr <- cov2cor(Sigmas.alr[[host]])
  # stack all unique pairwise correlations between microbes in a row associated with this host
  interaction_pairs.clr[i,] <- c(corr_mat.clr)
  interaction_pairs.alr[i,] <- c(corr_mat.alr)
  interaction_pairs.clr.tri[i,] <- corr_mat.clr[upper.tri(corr_mat.clr, diag=F)]
  interaction_pairs.alr.tri[i,] <- corr_mat.alr[upper.tri(corr_mat.alr, diag=F)]
}

# ======================================================================================================
#   plot heatmap over all components, blocked by family
# ======================================================================================================

# order of interaction_pairs.clr is (1x1) (1x3) ... (1xD) (2x1) (2x2) ...
# as pairs these are:                 1     2   ...   D   1*D+1 1*D+2 ...
# insert "spacers" for visualization
spacer_width <- 10
spaced.clr <- matrix(NA, nrow(interaction_pairs.clr), ncol(interaction_pairs.clr)+(models$D-1)*spacer_width)
for(i in 1:models$D) {
  offset.1_1 <- (i-1)*models$D + (i-1)*spacer_width + 1
  offset.2_1 <- offset.1_1 + models$D - 1
  offset.1_2 <- (i-1)*models$D + 1
  offset.2_2 <- offset.1_2 + models$D - 1
  spaced.clr[,offset.1_1:offset.2_1] <- interaction_pairs.clr[,offset.1_2:offset.2_2]
}
df <- gather_array(spaced.clr, "correlation", "host", "pair")
p <- ggplot(df, aes(pair, host)) +
  geom_tile(aes(fill = correlation)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred", na.value='#888888')
ggsave(file.path(relative_path,plot_dir,paste0("microbe_pair_correlations_",level,"_blocked.png")),
       p, units="in", dpi=150, height=5, width=20)

# collect unique interactions only for subsequent stuff
label_pairs.clr.tri <- label_pairs.clr[upper.tri(label_pairs.clr, diag=F)]
label_pairs.alr.tri <- label_pairs.alr[upper.tri(label_pairs.alr, diag=F)]

# ======================================================================================================
#   plot heatmap over just one component
#   how coherent is this family over all its interactions over all individuals?
# ======================================================================================================

if(level == "family") {
  idx <- 7 # family Muribaculaceae
  idx <- 22 # family Erysipelotrichaceae
  idx <- 31 # family Family XI_2II
  idx <- 35 # family Helicobacteraceae
  idx <- 40 # family Streptococcaceae
}
if(level == "genus") {
  idx <- 32 # genus Solobacterium
  idx <- 44 # genus Coriobacteriaceae UCG-003
  idx <- 69 # genus Lachnospiraceae ND3007 group
  idx <- 76 # genus [Eubacterium] oxidoreducens group
  idx <- 84 # genus Butyricicoccus
}

# plot heatmap
df <- gather_array(interaction_pairs.clr[,((idx-1)*models$D+1):(idx*models$D)], "correlation", "host", "pair")
p <- ggplot(df, aes(pair, host)) +
  geom_tile(aes(fill = correlation), colour = "white") +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave(file.path(relative_path,plot_dir,paste0("microbe_pair_correlations_",level,"_",idx,".png")),
       p, units="in", dpi=150, height=5, width=10)

# delete this once we're sure we don't need it

#retain_idx <- c()
#for(i in 1:ncol(interaction_pairs.clr.tri)) {
#  microbe_pair <- as.numeric(strsplit(label_pairs.clr.tri[i], "_")[[1]])
#  if(microbe_pair[1] == idx | microbe_pair[2] == idx) {
#    #cat("retaining:",microbe_pair[1],"x",microbe_pair[2],"\n")
#    retain_idx <- c(retain_idx, i)
#  }
#}
# subset to the interactions that include this component (should be D-1)
#subset_interaction_pairs.clr <- interaction_pairs.clr[,retain_idx]

# get distance of each pair
#d <- dist(t(subset_interaction_pairs.clr))
#clustering <- hclust(d)

# plot heatmap
#df <- gather_array(subset_interaction_pairs.clr, "correlation", "host", "pair")
#p <- ggplot(df, aes(pair, host)) +
#  geom_tile(aes(fill = correlation), colour = "white") +
#  scale_fill_gradient2(low = "darkblue", high = "darkred")
#ggsave(file.path(relative_path,plot_dir,paste0("microbe_pair_correlations_",level,"_",idx,".png")),
#       p, units="in", dpi=150, height=5, width=10)

# ======================================================================================================
#   gather some data and hierarchically cluster interactions
# ======================================================================================================

# average CLR abundance
counts <- otu_table(ref)@.Data
counts.clr <- clr(counts + 0.5)
avg_clr_abundance <- colMeans(counts.clr)

# get distance of each pair
d <- dist(t(interaction_pairs.clr.tri))
clustering <- hclust(d)

# reorder
interaction_pairs.clr.reordered <- interaction_pairs.clr.tri[,clustering$order]
label_pairs.clr.reordered <- label_pairs.clr.tri[clustering$order]

avg_interaction.clr <- colMeans(interaction_pairs.clr.reordered)


# here I'm using indices manually identified from
#   plot(avg_interaction.clr)
# to pull spans of strong negative, positive, and neural correlations
if(level == "genus") {
  neg_span <- 582:674
  pos_span <- 3193:3367
  neu_span <- 2500:2700
  whole_span <- 1:length(label_pairs.clr.reordered)
}

# ======================================================================================================
#   render NEGATIVE interaction biplots
# ======================================================================================================

if(level == "genus") {
  # get strong-ish average NEGATIVE interactions in the CLR
  interesting_idx.clr <- c(neg_span)
  df <- data.frame(x=c(), y=c(), sign=c(), value=c())
  for(p_idx in interesting_idx.clr) {
    interesting_pair <- label_pairs.clr.tri[p_idx]
    microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
    label.1 <- labels.clr[microbe_pair[1]]
    label.2 <- labels.clr[microbe_pair[2]]
    df <- rbind(df, data.frame(x=label.1, y=label.2, sign=sign(avg_interaction.clr[p_idx]),
                               value=2*abs(avg_interaction.clr[p_idx])))
    cat(paste0("Interesting pair (",p_idx,"): ",label.1,", ",label.2,"\n"))
  }
  
  text_df.left <- data.frame(y=1:length(unique(df$x)), label=unique(df$x))
  text_df.right <- data.frame(y=seq(1, length(unique(df$x)), length.out=length(unique(df$y))), label=unique(df$y))
  
  df <- cbind(df, node1=text_df.left$y[as.numeric(df$x)])
  df <- cbind(df, node2=text_df.right$y[as.numeric(df$y)])
  df$sign <- as.factor(df$sign)
  
  p <- ggplot() +
    geom_segment(data=df, aes(x=1, y=node1, xend=4, yend=node2, color=sign, size=value), alpha=0.2) +
    scale_color_manual(values=c("#E6194B")) +
    scale_size_identity() + # use the width specified by `weight`
    geom_text(data=text_df.left, aes(x=1, y=y, label=label), size=4) +
    geom_text(data=text_df.right, aes(x=4, y=y, label=label), size=4) +
    xlim(0, 5) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank()) 
  ggsave(file.path(relative_path,plot_dir,paste0("negative_interactions_",level,".png")),
         p, units="in", dpi=150, height=10, width=10)
}

# ======================================================================================================
#   render POSITIVE interaction biplots
# ======================================================================================================

if(level == "genus") {
  # get strong-ish average POSITIVE interactions in the CLR
  interesting_idx.clr <- c(pos_span)
  df <- data.frame(x=c(), y=c(), sign=c(), value=c())
  for(p_idx in interesting_idx.clr) {
    interesting_pair <- label_pairs.clr.tri[p_idx]
    microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
    label.1 <- labels.clr[microbe_pair[1]]
    label.2 <- labels.clr[microbe_pair[2]]
    df <- rbind(df, data.frame(x=label.1, y=label.2, sign=sign(avg_interaction.clr[p_idx]),
                               value=2*abs(avg_interaction.clr[p_idx])))
    cat(paste0("Interesting pair (",p_idx,"): ",label.1,", ",label.2,"\n"))
  }
  
  text_df.left <- data.frame(y=1:length(unique(df$x)), label=unique(df$x))
  text_df.right <- data.frame(y=seq(1, length(unique(df$x)), length.out=length(unique(df$y))), label=unique(df$y))
  
  df <- cbind(df, node1=text_df.left$y[as.numeric(df$x)])
  df <- cbind(df, node2=text_df.right$y[as.numeric(df$y)])
  df$sign <- as.factor(df$sign)
  
  p <- ggplot() +
    geom_segment(data=df, aes(x=1, y=node1, xend=4, yend=node2, color=sign, size=value), alpha=0.2) +
    scale_color_manual(values=c("#3CB44B")) +
    scale_size_identity() + # use the width specified by `weight`
    geom_text(data=text_df.left, aes(x=1, y=y, label=label), size=4) +
    geom_text(data=text_df.right, aes(x=4, y=y, label=label), size=4) +
    xlim(0, 5) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank()) 
  ggsave(file.path(relative_path,plot_dir,paste0("positive_interactions_",level,".png")),
         p, units="in", dpi=150, height=12, width=10)
}

# ======================================================================================================
#   plot heatmap over all components, hierarchically clustered
# ======================================================================================================

df <- gather_array(interaction_pairs.clr.reordered, "correlation", "host", "pair")
p <- ggplot(df, aes(pair, host)) +
  geom_tile(aes(fill = correlation)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave(file.path(relative_path,plot_dir,paste0("microbe_pair_correlations_",level,"_clustered.png")),
       p, units="in", dpi=150, height=5, width=15)

# make average relative abundance visualization
if(level == "genus") {
  neg_pairs <- label_pairs.clr.reordered[neg_span]
  neg_clr.1 <- c(); neg_clr.2 <- c()
  for(pair in neg_pairs) {
    microbe_pair <- as.numeric(strsplit(pair, "_")[[1]])
    temp.1 <- avg_clr_abundance[microbe_pair[1]]
    temp.2 <- avg_clr_abundance[microbe_pair[2]]
    if(temp.1 < temp.2) {
      neg_clr.1 <- c(neg_clr.1, temp.1)
      neg_clr.2 <- c(neg_clr.2, temp.2)
    } else {
      neg_clr.1 <- c(neg_clr.1, temp.2)
      neg_clr.2 <- c(neg_clr.2, temp.1)
    }
  }
  p <- ggplot(data.frame(x=1:length(neg_span), y=1, value=neg_clr.1), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_neg_1_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  p <- ggplot(data.frame(x=1:length(neg_span), y=1, value=neg_clr.2), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_neg_2_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  
  pos_pairs <- label_pairs.clr.reordered[pos_span]
  pos_clr.1 <- c(); pos_clr.2 <- c()
  for(pair in pos_pairs) {
    microbe_pair <- as.numeric(strsplit(pair, "_")[[1]])
    temp.1 <- avg_clr_abundance[microbe_pair[1]]
    temp.2 <- avg_clr_abundance[microbe_pair[2]]
    if(temp.1 < temp.2) {
      pos_clr.1 <- c(pos_clr.1, temp.1)
      pos_clr.2 <- c(pos_clr.2, temp.2)
    } else {
      pos_clr.1 <- c(pos_clr.1, temp.2)
      pos_clr.2 <- c(pos_clr.2, temp.1)
    }
  }
  p <- ggplot(data.frame(x=1:length(pos_span), y=1, value=pos_clr.1), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_pos_1_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  p <- ggplot(data.frame(x=1:length(pos_span), y=1, value=pos_clr.2), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_pos_2_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  
  neu_pairs <- label_pairs.clr.reordered[neu_span]
  neu_clr.1 <- c(); neu_clr.2 <- c()
  for(pair in neu_pairs) {
    microbe_pair <- as.numeric(strsplit(pair, "_")[[1]])
    temp.1 <- avg_clr_abundance[microbe_pair[1]]
    temp.2 <- avg_clr_abundance[microbe_pair[2]]
    if(temp.1 < temp.2) {
      neu_clr.1 <- c(neu_clr.1, temp.1)
      neu_clr.2 <- c(neu_clr.2, temp.2)
    } else {
      neu_clr.1 <- c(neu_clr.1, temp.2)
      neu_clr.2 <- c(neu_clr.2, temp.1)
    }
  }
  p <- ggplot(data.frame(x=1:length(neu_span), y=1, value=neu_clr.1), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_neu_1_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  p <- ggplot(data.frame(x=1:length(neu_span), y=1, value=neu_clr.2), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_neu_2_",level,".png")),
         p, units="in", dpi=150, height=1, width=5)
  
  whole_pairs <- label_pairs.clr.reordered[whole_span]
  whole_clr.1 <- c(); whole_clr.2 <- c()
  for(pair in whole_pairs) {
    microbe_pair <- as.numeric(strsplit(pair, "_")[[1]])
    temp.1 <- avg_clr_abundance[microbe_pair[1]]
    temp.2 <- avg_clr_abundance[microbe_pair[2]]
    if(temp.1 < temp.2) {
      whole_clr.1 <- c(whole_clr.1, temp.1)
      whole_clr.2 <- c(whole_clr.2, temp.2)
    } else {
      whole_clr.1 <- c(whole_clr.1, temp.2)
      whole_clr.2 <- c(whole_clr.2, temp.1)
    }
  }
  p <- ggplot(data.frame(x=1:length(whole_span), y=1, value=whole_clr.1), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_whole_1_",level,".png")),
         p, units="in", dpi=150, height=3, width=15)
  p <- ggplot(data.frame(x=1:length(whole_span), y=1, value=whole_clr.2), aes(x, y)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low="darkgreen", high="darkorange", limits=c(min(avg_clr_abundance), max(avg_clr_abundance)))
  ggsave(file.path(relative_path,plot_dir,paste0("colorbar_whole_2_",level,".png")),
         p, units="in", dpi=150, height=3, width=15)
}


# ======================================================================================================
#   sanity check strong interactions by plotting some strongly correlated ALR features
# ======================================================================================================

# get strong average interactions in the ALR
# plot the pairs across random individuals as a sanity check
avg_interaction.alr <- colMeans(interaction_pairs.alr)

interesting_idx.alr <- which(abs(avg_interaction.alr) >= 0.7)
for(p_idx in interesting_idx.alr) {
  interesting_pair <- label_pairs.alr.tri[p_idx]
  microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
  cat("Interesting pair (",p_idx,"): (",microbe_pair[1],")",labels.alr[microbe_pair[1]],", (",microbe_pair[2],")",labels.alr[microbe_pair[2]],"\n")
  sampled_hosts <- models$individuals[sample(1:length(models$individuals))[1:2]]
  for(host in sampled_hosts) {
    # sanity check universal pairs via predictive plots
    #plot_ribbons_individuals(c(host), level, timecourse=FALSE, covcor=FALSE, predict_coords=microbe_pair)
  }
}

# ======================================================================================================
#   from MAP estimates only (right now), look for microbes with a reasonably coherent set of
#   interactions across a large subset of individuals
# ======================================================================================================

# let's look at features matched at the same level of resolution; i.e. I'm not sure it makes sense
# to compare a microbial family to a group of other taxa resolved only to order
#resolved_idx <- as.vector(!is.na(tax[,level_no]))

host_prop_thresholds <- seq(0.5, 0.7, by=0.05)
interaction_prop_thresholds <- seq(0.5, 0.7, by=0.05)

if(level == "family") {
  host_prop_thresholds <- c(0.65)
  interaction_prop_thresholds <- c(0.65)
}
if(level == "genus") {
  host_prop_thresholds <- c(0.65)
  interaction_prop_thresholds <- c(0.6)
}

for(host_prop_threshold in host_prop_thresholds) {
  for(interaction_prop_threshold in interaction_prop_thresholds) {
    count <- 0
    for(i in 1:models$D) {
      offset.1 <- (i-1)*models$D + 1
      offset.2 <- offset.1 + models$D - 1
      temp <- interaction_pairs.clr[,offset.1:offset.2]
      temp2 <- apply(temp, c(1,2), sign)
      prop <- abs(colSums(temp2))/nrow(temp2)
      host_agreement <- sum(prop >= host_prop_threshold)/length(prop)
      if(host_agreement >= interaction_prop_threshold) {
        cat(paste0(labels.clr[i]," (",i,") is universal!\n"))
        count <- count + 1
      }
    }
    cat(host_prop_threshold,"x",interaction_prop_threshold,", count:",count,"\n")
  }
}







