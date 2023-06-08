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

evaluation_new <- function(combination, len_1, len_2){
  num_rows <- nrow(combination)
  num_cols <- ncol(combination)
  result_matrix <- matrix(0, nrow = num_rows, ncol = num_cols-1)
  
  
  # Iterate over each row of the matrix
  for (i in 1:num_rows) {
    row <- combination[i, ]
    row <- lapply(row, as.character)
    clock = 0
    
    print(paste("i is", i))
    for (j in 1:(num_cols-1)) {
      
      if (max(result_matrix[i, ]) <= threshold){
        clock = clock + 1
        result_matrix[i, clock] <- calculate_dimer(row[[j]], row[[num_cols]])$temp
      }
    }
  }
  
  result_matrix <- as.data.frame(result_matrix)
  row_indices <- which(apply(result_matrix, 1, function(row) all(row < threshold)))
  #print(row_indices)
  indices_1 <- unique(ceiling(row_indices/len_1))
  if (len_1 == 1){
    indices_1 <- 1
  }
  indices_2 <- unique(sapply(row_indices,
                             function(x) if(x%%len_2 == 0){return(len_2)}
                             else
                             {return(x%%len_2)} ))
  #print(indices_1)
  #print(indices_2)
  good_slection = c(list(indices_1), list(indices_2))
  return(good_slection)
}



get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}