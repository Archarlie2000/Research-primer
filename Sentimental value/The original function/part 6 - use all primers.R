library(tidyverse)
library(TmCalculator)
library(openPrimeR)
library(BiocManager)
library(plyr)
library(qpcR)
library(ggplot2)
library(stringi)

df <- read.csv("geenralprimerset.csv")

df <- df[c(2,3,4)]

####### creat mismatches

get_strong1 <- function(x){
  target <- str_sub(x , - 3, - 3) 
  target <- complement(target)
  
  if (target == "A") {temp <- "G"} else
    if (target == "G") {temp <- "A"} else
      if (target == "C") {temp <- "T"} else
        if (target == "T") {temp <- "C"}
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}

## Mismatching on Ts
get_strong2 <- function(x){
  target <- str_sub(x , - 3, - 3) 
  target <- complement(target)
  
  if (target == "T") {
    temp <- "T"
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)}
  else
    return(NULL)
}


get_medium1 <- function(x){
  target <- str_sub(x , - 3, - 3) 
  target <- complement(target)
  
  if (target == "A") {temp <- "A"} else
    if (target == "G") {temp <- "G"} else
      if (target == "C") {temp <- "C"} else
        return(NULL)
  
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}

get_weak1 <- function(x){
  target <- str_sub(x , - 3, - 3) 
  target <- complement(target)
  
  if (target == "C") {temp <- "A"} else
    if (target == "A") {temp <- "C"} else
      if (target == "G") {temp <- "T"} else
        if (target == "T") {temp <- "G"}
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}

# populate all the possible primers
mismatch_list <- df %>% 
  mutate(strong_mismatch_1 = map(forward, get_strong1),
         strong_mismatch_2 = map(forward, get_strong2),
         Medium_mismatch = map(forward, get_medium1),
         Weak_mismatch = map(forward, get_weak1))


# populate data frame for primer 3 for cross checking

Primer3_output <- mismatch_list[c(1,4,5,6,7,3)]

title <- rbind(cbind(paste(Primer3_output$name, "-","strong_mismatch_1")),
               cbind(paste(Primer3_output$name, "-","strong_mismatch_2")),
               cbind(paste(Primer3_output$name, "-","Medium_mismatch")),
               cbind(paste(Primer3_output$name, "-","Weak_mismatch")))


primers <- rbind(cbind(Primer3_output$strong_mismatch_1),
                 cbind(Primer3_output$strong_mismatch_2),
                 cbind(Primer3_output$Medium_mismatch),
                 cbind(Primer3_output$Weak_mismatch))

reverse <- rbind(cbind(mismatch_list$reverse),
                 cbind(mismatch_list$reverse),
                 cbind(mismatch_list$reverse),
                 cbind(mismatch_list$reverse))

Primer3_output <- cbind(title, primers, reverse)
colnames(Primer3_output) <- c("name", "forward", "reverse")
Primer3_output <- as.data.frame(Primer3_output) %>% 
  drop_na()


##### get Tm


Primer3_output <- cbind(Primer3_output, "", "", "", "")

for (i in 1:nrow(Primer3_output)){
    a <- Tm_NN(as.character(Primer3_output[i,3]), Na = 50, saltcorr = "SantaLucia1998-1")
    b <- Tm_NN(as.character(Primer3_output[i,2]), Na = 50, saltcorr = "SantaLucia1998-1")
    Primer3_output[i, 5] <- round(a$Tm, 2)
    Primer3_output[i, 4] <- round(b$Tm, 2)
}

colnames(Primer3_output) <- c("name", "forward", "reverse", "Tm1", "Tm2", "diff", "GC2")
Primer3_output <- as.data.frame(Primer3_output)


Primer3_output$diff <- abs(as.numeric(Primer3_output$Tm1) - 
                             as.numeric(Primer3_output$Tm2))
Primer3_output$GC2 <- (str_count(Primer3_output$reverse, "G") + 
                         str_count(Primer3_output$reverse, "C") 
              / str_length(Primer3_output$reverse) * 100) 

Primer3_output <- filter(Primer3_output,  diff < 2)

######## Use the original primer set to check individual compatibility


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



# Rename the column
outputframe <- data.frame(matrix(ncol = 6, nrow = 67*67))
colnames(outputframe) <- c("name", "forward", "reverse", "revcom", "template", "name2")



# Match every forward with every reverse primer for cross checking
k <- 0
for (i in 1:nrow(Primer3_output)){
  for (j in 1:nrow(Primer3_output)){
    k <- j+(i-1)*nrow(Primer3_output)
    outputframe[k,1] <- Primer3_output[i,1]
    outputframe[k,2] <- Primer3_output[i,2]
    outputframe[k,3] <- Primer3_output[j,3]
    outputframe[k,6] <- Primer3_output[j,1]
  }
}


outputframe <- read.csv("temp05")


#Extra sequncy
dymmy <- "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"


# Put name, template, reverse of complementary of reverse together
for (i in 1:nrow(outputframe)){
  print(i)
  outputframe[i,1] <- paste(outputframe[i,1], outputframe[i,6], sep  = "-")
  outputframe[i,4] <- complement(outputframe[i,3])
  outputframe[i,4] <- reverse_chars(outputframe[i,4])
  outputframe[i,5] <- gsub(" ", "", paste(outputframe[i,2], dymmy,
                                          outputframe[i,4], sep = ""))
}


#Enter output data frame into an input file for primer3
firstLine <- paste("SEQUENCE_ID=", outputframe$name, sep = "")
secondLine <- paste("SEQUENCE_TEMPLATE=", outputframe$template, sep = "")
thirdLine <- paste("SEQUENCE_PRIMER=", outputframe$forward, sep  = "")
fourthLine <- paste("SEQUENCE_PRIMER_REVCOMP=", outputframe$reverse, sep  = "")
otherLines <- "PRIMER_TASK=generic
PRIMER_OPT_SIZE=21
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=30
PRIMER_PRODUCT_SIZE_RANGE=30-500
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
="

finalOutput <- paste(firstLine, secondLine, thirdLine, fourthLine, otherLines, sep = "\n")

write(finalOutput, "input310.txt")



##### get data for mismatching set





