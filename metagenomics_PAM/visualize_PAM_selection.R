# this files visualizes selected PAM samples as a sanity check
# we want these well distributed among the full sample set with respect to:
#   (1) time
#   (2) dynamics distance (at least in terms of filtered individuals; this is all we can get)
#   (3) phylogenetic distance (e.g. weighted UniFrac)

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

kk <- c(500, 1000, 1500, 2000)

for(ll in c("genus", "family", "phylum")) {
  # load full data
  data <- load_and_filter(level=ll)
  metadata <- sample_data(data)

  for(k in kk) {
    cat("Level:",ll,", sample number:",k,"\n")
    samples <- read.table(file=file.path(relative_path,output_dir,paste0('PAM_selected_samples_',k,'.tsv')), sep='\t', header=TRUE, row.names=1)

    # `samples` has the format:
    #                        SID collection_date sname sample_status
    #   1 11408-ACGCGTGCCCGG-395      2001-01-02   DUI             0
    #   2 12053-GAACCAGTACTC-409      2001-01-03   ACA             0
    #   3 11412-CATTCGTGGCGT-397      2001-01-25   DUI             0
    #   4 11412-GTGCGGTTCACT-397      2001-01-29   VEX             0
    #   5 12053-GAAATCTTGAAG-409      2001-01-29   FAX             0
    #   6 12053-CTCGCCCTCGCC-409      2001-02-03   NAP             0

    # visualize w/r/t TIME

    if(ll == "genus") { # just do this once
      collection_date <- metadata$collection_date
      sname <- metadata$sname
      collection_date <- as.POSIXct(collection_date)
      df <- data.frame(date=collection_date, sname=as.factor(sname), selected=FALSE)

      df <- rbind(df, data.frame(date=as.POSIXct(samples$collection_date), sname=samples$sname, selected=TRUE))
      p <- ggplot(df) +
        geom_point(aes(x=date, y=sname, color=selected)) +
        scale_x_datetime(date_breaks="1 month", labels=date_format("%Y-%m-%d")) +
        theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12,face="bold")) +
        scale_colour_manual(values = c("gray", "orange"))
      ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_time_",k,".png")), p, units="in", dpi=200, height=10, width=20)
    }

    # visualize w/r/t dynamics distance

    # to approximate this we will:
    #   (1) load distances between individual MAP estimates in terms of dynamics
    #   (2) indicate selected individuals by color
    #   (3) indicate number of samples chosen from each individual by size of geom_point() dot
    indiv_dist <- readRDS(file.path(relative_path,output_dir,paste0("saved_distance_Sigma_",ll,"_MAP.rds")))
    embedding <- cmdscale(indiv_dist, k=20)
    df <- data.frame(sname=sname_list, x=embedding[,1], y=embedding[,2],
                                       x2=embedding[,9], y2=embedding[,10],
                                       x3=embedding[,19], y3=embedding[,20], selected=FALSE)
    for(i in 1:nrow(samples)) {
      idx <- which(df$sname == as.character(samples$sname[i]))
      df[idx,]$selected <- TRUE
    }
    p <- ggplot(df) +
      geom_point(aes(x=x, y=y, color=selected), size=2) +
      scale_colour_manual(values = c("gray", "orange")) +
      xlab("PCoA 1 (Riemannian distance)") +
      ylab("PCoA 2 (Riemannian distance)")
    ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_dynamics_",k,"_1x2.png")), p, units="in", dpi=150, height=10, width=10)
    p <- ggplot(df) +
      geom_point(aes(x=x2, y=y2, color=selected), size=2) +
      scale_colour_manual(values = c("gray", "orange")) +
      xlab("PCoA 9 (Riemannian distance)") +
      ylab("PCoA 10 (Riemannian distance)")
    ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_dynamics_",k,"_9x10.png")), p, units="in", dpi=150, height=10, width=10)
    p <- ggplot(df) +
      geom_point(aes(x=x3, y=y3, color=selected), size=2) +
      scale_colour_manual(values = c("gray", "orange")) +
      xlab("PCoA 19 (Riemannian distance)") +
      ylab("PCoA 20 (Riemannian distance)")
    ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_dynamics_",k,"_19x20.png")), p, units="in", dpi=150, height=10, width=10)

    if(ll == "family") { # just do this once

      # visualize w/r/t whole data set in Aitchison distance
      
      SIDs_all <- rownames(metadata)
      selected_idx <- SIDs_all %in% samples$SID

      cat("Calculating Aitchison distance at level:",ll,"\n")
      clr_all <- driver::clr(otu_table(data)@.Data + 0.5) # expects samples as rows; returns samples as rows
      clr_dist <- dist(clr_all)
      clr_embed <- cmdscale(clr_dist)
      df <- data.frame(sname=SIDs_all, x=clr_embed[,1], y=clr_embed[,2], selected=FALSE)
      df[selected_idx,]$selected <- TRUE

      p <- ggplot(df) +
        geom_point(aes(x=x, y=y, color=selected), size=2) +
        scale_colour_manual(values = c("gray", "orange")) +
        xlab("PCoA 1 (Aitchison distance)") +
        ylab("PCoA 2 (Aitchison distance)")
      ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_Aitchison_",k,".png")), p, units="in", dpi=150, height=10, width=10)

      # visualize w/r/t phylogenetic distance (weighted UniFrac)

      cat("Calculating UniFrac distance at level:",ll,"\n")

      # for reasonable running time, try this will 1/4 of the samples
      retain_prop <- 0.1

      if(file.exists(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_selected_sids.rds")) &
        file.exists(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_all_sids.rds"))) {
        sub_sid_selected <- readRDS(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_selected_sids.rds"))
        sub_sid_all <- readRDS(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_all_sids.rds"))
      } else {
        retain_from_selected <- as.logical(rbinom(k, 1, retain_prop))
        retain_from_all <- as.logical(rbinom(length(selected_idx), 1, retain_prop))
        sub_sid_selected <- as.character(samples$SID[retain_from_selected])
        sub_sid_all <- sort(unique(c(sub_sid_selected, SIDs_all[retain_from_all])))
        saveRDS(sub_sid_selected, paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_selected_sids.rds"))
        saveRDS(sub_sid_all, paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_all_sids.rds"))
      }

      sub_data <- subset_samples(data, sample_id %in% sub_sid_all)
      
      cat("\tCalculating UniFrac distance over",length(sub_sid_all),"samples...\n")
      if(file.exists(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_dist.rds"))) {
        unifrac_dist <- readRDS(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_dist.rds"))
      } else {
        unifrac_dist <- UniFrac(sub_data, weighted=TRUE, parallel=TRUE)
      }
      cat("\tEmbedding UniFrac distances over",length(sub_sid_all),"samples...\n")
      if(file.exists(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_embed.rds"))) {
        unifrac_embed <- readRDS(paste0("output/PAM_UniFrac_",k,"_",retain_prop,"_embed.rds"))
      } else {
        unifrac_embed <- cmdscale(unifrac_dist)
      }

      df <- data.frame(sname=sub_sid_all, x=unifrac_embed[,1], y=unifrac_embed[,2], selected=FALSE)
      df[df$sname %in% sub_sid_selected,]$selected <- TRUE

      # all samples; this takes an eternity
      # unifrac_dist <- UniFrac(data, weighted=TRUE, parallel=TRUE)
      # unifrac_embed <- cmdscale(unifrac_dist)
      # df <- data.frame(sname=SIDs_all, x=unifrac_embed[,1], y=unifrac_embed[,2], selected=FALSE)
      # df[selected_idx,]$selected <- TRUE

      p <- ggplot(df) +
        geom_point(aes(x=x, y=y, color=selected), size=2) +
        scale_colour_manual(values = c("gray", "orange")) +
        xlab("PCoA 1 (UniFrac distance)") +
        ylab("PCoA 2 (UniFrac distance)")
      ggsave(file.path(relative_path,output_dir,paste0("PAM_",ll,"_diagnostic_UniFrac_",k,"_",retain_prop,"_1x2.png")), p, units="in", dpi=150, height=10, width=10)

    }
  }
}

