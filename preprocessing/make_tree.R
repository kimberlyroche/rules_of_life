# this file is lifted from a tutorial online and builds a tree associated with the sequence content
# in the existing ABRP phyloseq files I've got; we'll resave the processed data files with a `_tree.rds`
# extension
#
# note: a tree is necessary if we want to (elsewhere) calculate UniFrac distance

relative_path <- ".."

source(file.path(relative_path,"include/R/general.R"))

library(dada2)
library(phangorn)
library(msa)
library(DECIPHER)

set.seed(4)

load(file.path(relative_path,data_dir,"glom_data_genus_reps.RData")) # loads glom_data
ftbl <- as(otu_table(glom_data), "matrix") # your feature table

seqs <- getSequences(ftbl)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit <- pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)

print(Sys.time())
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
print(Sys.time())
detach("package:phangorn", unload=TRUE)

print(is.rooted(phy_tree(fitGTR$tree)))

phy_tree(glom_data) <- phy_tree(fitGTR$tree)
saveRDS(glom_data, file.path(relative_path,data_dir,"glom_data_genus_reps_tree.rds"))

