library(tidyverse)
library(TmCalculator)
library(openPrimeR)
library(BiocManager)
library(plyr)
library(qpcR)
library(ggplot2)
library(stringi)



primers <- read.delim("output_2021.txt", sep= '=', header = FALSE)


primersFiltered <- filter(primers, V1 %in% c('SEQUENCE_ID', 
                                             'SEQUENCE_PRIMER',
                                             'SEQUENCE_PRIMER_REVCOMP',
                                             'PRIMER_LEFT_0_HAIRPIN_TH',
                                             'PRIMER_LEFT_0_END_STABILITY',
                                             'PRIMER_RIGHT_0_END_STABILITY'))


getname <- data.frame(name = primersFiltered$V2[primersFiltered$V1 == "SEQUENCE_ID"])
getseq2 <- data.frame(seq2 = primersFiltered$V2[primersFiltered$V1 == "SEQUENCE_PRIMER"])
getseq3 <- data.frame(seq3 = primersFiltered$V2[primersFiltered$V1 == "SEQUENCE_PRIMER_REVCOMP"])
gethair <- data.frame(seq3 = primersFiltered$V2[primersFiltered$V1 == "PRIMER_LEFT_0_HAIRPIN_TH"])
sta1 <- data.frame(seq3 = primersFiltered$V2[primersFiltered$V1 == "PRIMER_LEFT_0_END_STABILITY"])
sta2 <- data.frame(seq3 = primersFiltered$V2[primersFiltered$V1 == "PRIMER_RIGHT_0_END_STABILITY"])


getcross <- qpcR:::cbind.na(getname$name, getseq2$seq2, 
                            getseq3, gethair, sta1, sta2)  


colnames(getcross) <- c("name", "forward", "reverse", "hairpin", "stability_1", "stability_2")

nonreplica <- distinct(getcross, forward, reverse, .keep_all= TRUE)
nonreplica <- as.data.frame(nonreplica)


