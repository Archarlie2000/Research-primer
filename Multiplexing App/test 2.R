# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rprimer")
# 
# BiocManager::install("Biostrings")

library(rprimer)
library(Biostrings)

system.file("extdata", "example_alignment.txt", package = "rprimer")
filepath <- system.file("extdata", "example_alignment.txt", package = "rprimer")



myAlignment <- readDNAMultipleAlignment(filepath, format = "fasta")



myMaskedAlignment <- myAlignment
colmask(myMaskedAlignment, invert = TRUE) <- c(3000:4000, 5000:6000)

myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)

plotData(myConsensusProfile)

myOligos <- designOligos(myConsensusProfile)
