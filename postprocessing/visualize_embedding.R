# this file labels the already-calculated embedding by the specified annotation (e.g. host, group...)

relative_path <- ".."

source(file.path(relative_path,"include/R/GP.R"))

# available labels
#   individual
#   group
#   matgroup
#   counts
#   density
#   mom
#   dad
#   momrank
#   drought
#   largegroup
#   momdied
#   competingsib
#   earlyadversity
#   birthrate_all
#   birthrate_surviving
#   alphadiv
#   wetproportion
#   metagenomics -- note: this is a list of individuals selected for metagenomics sample collection

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) {
  stop("Arguments: (level) (measure=Sigma|Lambda) (label type) (opt: sname list)", call.=FALSE)
}
level <- args[1]
which_measure <- args[2]
label_type <- args[3]
indiv_list <- NULL
if(label_type == "metagenomics" & length(args) > 3) {
  indiv_list <- args[4:length(args)]
}

plot_ordination(level, which_measure, label_type, legend=TRUE, sname_list=indiv_list)
