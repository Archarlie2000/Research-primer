library(biomaRt)
library(tidyverse)
library(stringr)
library(sys)
library(TmCalculator)
library(stringi)


# list of snps we want information from
# would be good to read this in from a file in the future
snp_list <- c("rs25", "rs16944", "rs1884", "rs17287498")

# the length of flanking sequences we will retrieve
upStream <- c("500")
downStream <- c("500")

# use biomaRt to connect to the human snp database
snpmart <- useMart("ENSEMBL_MART_SNP",dataset = "hsapiens_snp")

# this actually goes and gets the flanking sequences for our snps in the snp_list
snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'), 
                      filters = c('snp_filter', 'upstream_flank', 'downstream_flank'), 
                      checkFilters = FALSE, 
                      values = list(snp_list, upStream, downStream), 
                      mart = snpmart, 
                      bmHeader = TRUE)


# for each sequence, there is a portion that lists the different bases possible
# the format of this portion of the sequence is %N/N%

# this splits up the sequence into the upstream flank, variants in the format N/N, downstream flank
snp_sequence_split <- as.data.frame(str_split(snp_sequence$`Variant sequences`, "%"))

# the data frame just created has each snp in a column, need to name the columns
colnames(snp_sequence_split) <- snp_sequence$`Variant name`
# name the rows of the data frame
rownames(snp_sequence_split) <- c("upstream", "variants", "downstream")

# want to transpose the data frame so each variant is a row instead of a column
snps <- t(snp_sequence_split)

# now work on getting the variants split
snpsTibble <- as_tibble(snps, rownames = NA)

vars_split <- str_split(snpsTibble$variants, "/")
# vars_split now contains the possible variants, but there are possibly different numbers of variants
# need to make all entries in vars_split the same length
#n <- max(lengths(vars_split))
# set n to be the most possible variants we could possibly see at any snp position
# we'll go with 10, which should be much larger than any number we'll see
# can get rid of extras later
n <- 10
vars_split_uniform <- as.data.frame(lapply(vars_split, `length<-`, n))

# vars_split_uniform now has each snpID in a column and the possible variants as the rows
# need to name things to keep that straight
colnames(vars_split_uniform) <- snp_sequence$`Variant name`

vars_split_transposed <- t(vars_split_uniform)
variantsTibble <- as_tibble(vars_split_transposed, rownames = NA)

varListFinal <- mutate(variantsTibble, snpsTibble)
varListFinal['snpID'] <- snp_sequence$`Variant name`

# make it a tibble
variantsTibbleFinal <- as_tibble(varListFinal, rownames = NA)

variantsTibbleFinal2 <- pivot_longer(variantsTibbleFinal, 
                                     cols = V1:V10,
                                     names_to = "variations", 
                                     values_to = "observations")

# we padded the length of the variants columns earlier to make sure we 
# could handle any length we might see in the pivot_longer
# now, simply remove any rows that contain NA in the observations column
variantsTrimmed <- drop_na(variantsTibbleFinal2)

# variantsTrimmed has everything we need to make the output for calling primer3
# need to reorder then combine some columns
# we want one string that is the upstream sequence, then the observation of the snp
# then the downstream sequence all together in one string
variantsTrimmed <- variantsTrimmed %>% relocate(snpID)
variantsTrimmed <- variantsTrimmed %>% relocate(observations, .before = variants)
variantsTrimmed <- variantsTrimmed %>% relocate(downstream, .before = variants)
variantsTrimmed <- variantsTrimmed %>% unite("sequence", upstream:downstream, sep = "")


# add columns for the substrings leading up to and including the variant site
variantsTrimmed <- variantsTrimmed %>% mutate(left18 = str_sub(sequence, 484, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left19 = str_sub(sequence, 483, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left20 = str_sub(sequence, 482, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left21 = str_sub(sequence, 481, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left22 = str_sub(sequence, 480, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left23 = str_sub(sequence, 479, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left24 = str_sub(sequence, 478, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left25 = str_sub(sequence, 477, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left26 = str_sub(sequence, 476, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left27 = str_sub(sequence, 475, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left28 = str_sub(sequence, 474, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left29 = str_sub(sequence, 473, 501))
variantsTrimmed <- variantsTrimmed %>% mutate(left30 = str_sub(sequence, 472, 501))


# pivot longer so each left primer gets on it's own row in the tibble
variantsTrimmed2 <- pivot_longer(variantsTrimmed, 
                                     cols = left18:left30,
                                     values_to = "leftPrimers")

# calculate the Tm for the primers we will be specifying

Tm_NN(testPrimer, Na = 50, saltcorr = "SantaLucia1998-1")


firstLine <- paste("SEQUENCE_ID=", variantsTrimmed2$snpID, variantsTrimmed2$name, sep = "")
secondLine <- paste("SEQUENCE_TEMPLATE=", variantsTrimmed2$sequence, sep = "")
thirdLine <- paste("SEQUENCE_PRIMER=", variantsTrimmed2$leftPrimers, sep  = "")

otherLines <- "PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=21
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=30
PRIMER_PRODUCT_SIZE_RANGE=50-500
PRIMER_EXPLAIN_FLAG=1
="

finalOutput <- paste(firstLine, secondLine, thirdLine, otherLines, sep = "\n")

write(finalOutput, "input.txt")

exec_wait("pwd")
exec_wait("ls")
outputPrimer3 <- exec_wait("primer3-main/src/primer3_core", args = "input.txt", std_out = "output22.txt")

