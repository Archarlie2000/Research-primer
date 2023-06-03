get_strong1 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "A") {temp <- "G"} else
    if (target == "G") {temp <- "A"} else
      if (target == "C") {temp <- "T"} else
        if (target == "T") {temp <- "C"}
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}

## Get strong mismatch for the first three bp
left_flanking_get_strong1 <- function(x){
  temp <- ""
  target <- str_sub(x , 3, 3)
  
  if (target == "A") {temp <- "G"} else
    if (target == "G") {temp <- "A"} else
      if (target == "C") {temp <- "T"} else
        if (target == "T") {temp <- "C"}
  substring(x, 3, 3) <- temp
  return(x)
}

## Get strong mismatch for the last three bp, there are two types
get_strong2 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "T") {
    temp <- "T"
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)}
  else
    return("N")
}

## Get strong mismatch for the first three bp
left_flanking_get_strong2 <- function(x){
  temp <- ""
  target <- str_sub(x , 3, 3)
  if (target == "T") {
    temp <- "T"
    substring(x, 3, 3) <- temp
    return(x)}
  else
    return("N")
}


## Get medium mismatch for the last three bp
get_medium1 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "A") {temp <- "A"} else
    if (target == "G") {temp <- "G"} else
      if (target == "C") {temp <- "C"} else
        return("N")
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
  
}

## Get weak mismatch for the first three bp
left_flanking_get_medium1 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "A") {temp <- "A"} else
    if (target == "G") {temp <- "G"} else
      if (target == "C") {temp <- "C"} else
        return("N")
  substring(x, 3, 3) <- temp
  return(x)
}


## Get weak mismatch for the last three bp
get_weak1 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "C") {temp <- "A"} else
    if (target == "A") {temp <- "C"} else
      if (target == "G") {temp <- "T"} else
        if (target == "T") {temp <- "G"}
  substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
}


## Get weak mismatch for the first three bp
left_flanking_get_weak1 <- function(x){
  temp <- ""
  target <- str_sub(x , - 3, - 3)
  
  if (target == "C") {temp <- "A"} else
    if (target == "A") {temp <- "C"} else
      if (target == "G") {temp <- "T"} else
        if (target == "T") {temp <- "G"}
  substring(x, 3, 3) <- temp
  return(x)
}


## Get reversed sequnce of a string (There is a beter way to do this, but)
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


## Warngling SNP list to individual rows
list_seq <- function(snp) {
  first_position <- unlist(gregexpr('/', snp))[1]
  number_slash <- str_count(snp, "/")
  block <- str_sub(snp, first_position -1 , 
                   first_position - 2 + number_slash*2 + 1) 
  block <- gsub("/", "", block)
  block <- strsplit(block, "")
  
  for (i in block) {
    k <- paste(str_sub(snp, 1, first_position - 2),
               i,
               str_sub(snp, first_position - 2 + number_slash * 2 + 2, str_length(snp) ),
               sep = "")
    
  }
  k = gsub("%", "", k)
  k = k[!grepl("W", k)]
  return(k)
}