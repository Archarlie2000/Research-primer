library(rprimer)
library(Biostrings)

system.file("extdata", "example_alignment.txt", package = "rprimer")
filepath <- system.file("extdata", "example_alignment.txt", package = "rprimer")
myAlignment <- readDNAMultipleAlignment(filepath, format = "fasta")

## Mask everything but position 3000 to 4000 and 5000 to 6000
myMaskedAlignment <- myAlignment

colmask(myMaskedAlignment, invert = TRUE) <- c(3000:4000, 5000:6000)

myConsensusProfile <- consensusProfile(myAlignment, ambiguityThreshold = 0.05)
plotData(myConsensusProfile)

selection <- myConsensusProfile[
  myConsensusProfile$position >= 5000 & myConsensusProfile$position <= 5800, 
]

plotData(selection)

myOligos <- designOligos(myConsensusProfile)
