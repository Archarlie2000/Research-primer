---
title: "CASE STUDY TITLE"
author: "YOUR NAME"
date: "April 08, 2023"
output:
  html_document:  
    keep_md: true
    toc: true
    toc_float: true
    code_folding: hide
    fig_height: 6
    fig_width: 12
    fig_align: 'center'
---




```r
library(reticulate)
library(biomaRt)
library(tidyverse)
library(stringr)
library(sys)
library(TmCalculator)
library(stringi)
```



```r
snp_list <- c("rs25", "rs16944", "rs1884", "rs17287498")

upStream <- c("500")
downStream <- c("500")

snpmart <- useMart("ENSEMBL_MART_SNP",dataset = "hsapiens_snp")

# this actually goes and gets the flanking sequences for our snps in the snp_list
snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'), 
                      filters = c('snp_filter', 'upstream_flank', 'downstream_flank'), 
                      checkFilters = FALSE, 
                      values = list(snp_list, upStream, downStream), 
                      mart = snpmart, 
                      bmHeader = TRUE)


snp_sequence_split <- as.data.frame(str_split(snp_sequence$`Variant sequences`, "%"))

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
variantsTrimmed <- variantsTrimmed %>% relocate(observations, .before = variations)
variantsTrimmed <- variantsTrimmed %>% relocate(downstream, .before = variations)
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
```






```python


Big_template = r.variantsTrimmed2


import numpy
import primer3 as bindings
import pandas as pd


print(Big_template)

```

```
##        snpID  ...                     leftPrimers
## 0       rs25  ...              ATCAAATTCCAATTGCAT
## 1       rs25  ...             TATCAAATTCCAATTGCAT
## 2       rs25  ...            CTATCAAATTCCAATTGCAT
## 3       rs25  ...           ACTATCAAATTCCAATTGCAT
## 4       rs25  ...          GACTATCAAATTCCAATTGCAT
## ..       ...  ...                             ...
## 138  rs16944  ...      TACCTTGGGTGCTGTTCTCTGCCTCA
## 139  rs16944  ...     CTACCTTGGGTGCTGTTCTCTGCCTCA
## 140  rs16944  ...    TCTACCTTGGGTGCTGTTCTCTGCCTCA
## 141  rs16944  ...   CTCTACCTTGGGTGCTGTTCTCTGCCTCA
## 142  rs16944  ...  TCTCTACCTTGGGTGCTGTTCTCTGCCTCA
## 
## [143 rows x 5 columns]
```





