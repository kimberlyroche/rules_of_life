library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(phangorn); packageVersion("phangorn")
library(msa); packageVersion("msa")
library(DECIPHER); packageVersion("DECIPHER")
set.seed(4)

load("data/glom_data_genus_reps.RData") # loads glom_data
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
saveRDS(glom_data, "data/glom_data_genus_reps_tree.rds") 

