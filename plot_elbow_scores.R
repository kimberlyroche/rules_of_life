library(ClusterR)
source("include/R/general.R") # output_dir

# if (1) clustering by PAM for each taxonomic level over many k and (2) clean up and deposition
# of those scores into a file `output/all_elbow_scores.txt` has been done, plot these scores

scores <- read.table(file='output/all_elbow_scores.txt', sep='\t', header=TRUE)

p <- ggplot(scores) +
      geom_point(aes(x=clusters, y=score, color=level), size=2) +
      ylim(0, max(scores$score))
      #theme(legend.position = "none")
ggsave(paste0(output_dir,"elbow_intracluster_scores.png"), units="in", dpi=100, height=10, width=10)


