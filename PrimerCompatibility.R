library(tidyverse)
library(TmCalculator)
library(openPrimeR)
library(BiocManager)
library(plyr)
library(qpcR)
library(ggplot2)
library(stringi)



# need to read in the primer3 output file with potential primers
primers <- read.delim("output.txt", sep= '=', header = FALSE)

# first may want to filter this big table to get the primers for each snp into different tables
# try filtering into a new dataframe for each snp in a list

# I only want the sequence of the left primer and the possible right primers
primersFiltered <- primers[primers$V1 %in% c('SEQUENCE_ID', 
                                             'SEQUENCE_PRIMER', 
                                             'PRIMER_RIGHT_0_SEQUENCE',
                                             'PRIMER_RIGHT_0_GC_PERCENT',
                                             'PRIMER_RIGHT_0_END_STABILITY',
                                             'PRIMER_RIGHT_0_TM')]


#Recontructing the dataframe for Tm extraction
getname <- data.frame(name = primersFiltered$V2[primersFiltered$V1 == "SEQUENCE_ID"])
getseq1 <- data.frame(seq1 = primersFiltered$V2[primersFiltered$V1 == "SEQUENCE_PRIMER"])
getseq2 <- data.frame(seq2 = primersFiltered$V2[primersFiltered$V1 == "PRIMER_RIGHT_0_SEQUENCE"])

getTM <- qpcR:::cbind.na(getname$name, getseq1$seq1, getseq2$seq2, 
                         c("Na"), c("Na"), c("Na"), c("Na"))  

for (i in 1:nrow(getTM)){
  if (!is.na(getTM[i,3])){
  a <- Tm_NN(as.character(getTM[i,3]), Na = 50, saltcorr = "SantaLucia1998-1")
  b <- Tm_NN(as.character(getTM[i,2]), Na = 50, saltcorr = "SantaLucia1998-1")
  getTM[i, 5] <- round(a$Tm, 2)
  getTM[i, 4] <- round(b$Tm, 2)
  
  }
}

colnames(getTM) <- c("name", "forward", "reverse", "Tm1", "Tm2", "diff", "GC2")
getTM <- as.data.frame(getTM)

#Calculate the difference between mT
getTM$diff <- abs(as.numeric(getTM$Tm1) - as.numeric(getTM$Tm2))
getTM$GC2 <- (str_count(getTM$reverse, "G") + str_count(getTM$reverse, "C") 
              / str_length(getTM$reverse) * 100) 



#Remove primers that do not have a compatible reverse primers
getTMFiltered <- na.omit(getTMFiltered)
getTMFiltered <- filter(getTM,  diff < 2)

# the primers data frame has a column for the attributes and a column for the values,
# with a row for each attribute
# There are lots of attributes but the only ones I want to keep are:
# SEQUENCE_ID
# SEQUENCE_PRIMER
# N in the following lines will be a number (default 0 through 4) 
# PRIMER_RIGHT_N_SEQUENCE
# PRIMER_RIGHT_N this one is the position within the template and length of the primer
# PRIMER_RIGHT_N_TM
# PRIMER_RIGHT_N_GC_PERCENT
# PRIMER_RIGHT_N_SELF_ANY_TH
# PRIMER_RIGHT_N_SELF_END_TH
# PRIMER_RIGHT_N_HAIRPIN_TH
# PRIMER_RIGHT_N_END_STABILITY


outputframe <- qpcR:::cbind.na(getTMFiltered$name, getTMFiltered$forward, getTMFiltered$reverse, 
                         c("Na"), c("Na"), c("Na"), c("Na"))  


colnames(outputframe) <- c("name", "forward", "reverse", "k1", "k2", "k3", "k4")

reverse_chars <- function(string)
{
  # split string by characters
  string_split = strsplit(string, split = "")
  # reverse order
  rev_order = nchar(string):1
  # reversed characters
  reversed_chars = string_split[[1]][rev_order]
  # collapse reversed characters
  paste(reversed_chars, collapse = "")
} 






outputframe <- data.frame(matrix(ncol = 5, nrow = 420))
colnames(outputframe) <- c("name", "forward", "reverse", "k1", "k2")




k <- 0
for (i in 1:nrow(getTMFiltered)){
  for (j in 1:nrow(getTMFiltered)){
    k <- j+(i-1)*nrow(getTMFiltered)
    outputframe[k,1] <- getTMFiltered[i,1]
    outputframe[k,2] <- getTMFiltered[i,2]
    outputframe[k,3] <- getTMFiltered[j,3]
  }
}



for (i in 1:nrow(outputframe)){
  print(i)
  outputframe[i,4] <- complement(outputframe[i,3])
  outputframe[i,4] <- reverse_chars(outputframe[i,4])
  outputframe[i,5] <- gsub(" ", "", paste(outputframe[i,2], outputframe[i,3], sep = ""))
}

firstLine <- paste("SEQUENCE_ID=", outputframe$name, sep = "")
secondLine <- paste("SEQUENCE_TEMPLATE=", outputframe$k2, sep = "")
thirdLine <- paste("SEQUENCE_PRIMER=", outputframe$forward, sep  = "")
fourthLine <- paste("SEQUENCE_PRIMER_REVCOMP=", outputframe$k1, sep  = "")
otherLines <- "PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=0
PRIMER_PICK_RIGHT_PRIMER=0
PRIMER_OPT_SIZE=21
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=30
PRIMER_PRODUCT_SIZE_RANGE=50-500
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
="

finalOutput <- paste(firstLine, secondLine, thirdLine, fourthLine, otherLines, sep = "\n")

write(finalOutput, "input22.txt")
