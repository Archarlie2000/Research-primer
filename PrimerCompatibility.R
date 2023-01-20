library(tidyverse)
library(TmCalculator)
library(openPrimeR)
library(rmelting)
library(BiocManager)

BiocManager::install("rmelting")

# need to read in the primer3 output file with potential primers
primers <- read.delim("output.txt", sep= '=', header = FALSE)

# first may want to filter this big table to get the primers for each snp into different tables
# try filtering into a new dataframe for each snp in a list

# I only want the sequence of the left primer and the possible right primers
primersFiltered <- primers[primers$V1 %in% c('SEQUENCE_ID', 'SEQUENCE_PRIMER', 'PRIMER_RIGHT_0_SEQUENCE'), ]


testframe = data.frame(cola = c(primersFiltered[6,2], primersFiltered[2,2]))

testPrimer <- primersFiltered[6,2]

Tm_NN(ntseq = testframe$cola, Na = 50, saltcorr = "SantaLucia1998-1")
Tm_NN(ntseq = primersFiltered[2,2], Na = 50, saltcorr = "SantaLucia1998-1")




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

