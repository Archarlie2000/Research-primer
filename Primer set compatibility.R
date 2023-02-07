library(tidyverse)
library(TmCalculator)
library(openPrimeR)
library(BiocManager)
library(plyr)
library(qpcR)
library(ggplot2)
library(stringi)


primers <- read.delim("output.txt", sep= '=', header = FALSE)


primersFiltered <- primers[primers$V1 %in% c('SEQUENCE_ID', 
                                             'SEQUENCE_PRIMER', 
                                             'PRIMER_RIGHT_0_SEQUENCE',
                                             'PRIMER_RIGHT_0_GC_PERCENT',
                                             'PRIMER_RIGHT_0_END_STABILITY',
                                             'PRIMER_RIGHT_0_TM')]