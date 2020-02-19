# this file plots the score from `cluster_elbow_plot.R`

relative_path <- ".."

library(ClusterR)

source(file.path(relative_path,"include/R/general.R"))

# note: file {output_dir}/all_elbow_scores.txt must exist!

scores <- read.table(file=file.path(relative_path,output_dir,"all_elbow_scores.txt"), sep='\t', header=TRUE)

p <- ggplot(scores) +
      geom_point(aes(x=clusters, y=score, color=level), size=2) +
      ylim(0, max(scores$score))
      #theme(legend.position = "none")
ggsave(file.path(relative_path,output_dir,"elbow_intracluster_scores.png"), units="in", dpi=100, height=10, width=10)

