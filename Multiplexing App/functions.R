get_strong1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)
  
  if (target == "A") {temp <- "G"} else
    if (target == "G") {temp <- "A"} else
      if (target == "C") {temp <- "T"} else
        if (target == "T") {temp <- "C"}
  
  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  
  return(x)
}

get_strong2 <- function(x, type){
  temp <- ""
  
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)
  
  if (target == "T") {
    temp <- "T"
    if(type){
      substring(x, 3, 3) <- temp
    }else
      substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)}
  else
    return("N")
}

get_medium1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)
  
  if (target == "A") {temp <- "A"} else
    if (target == "G") {temp <- "G"} else
      if (target == "C") {temp <- "C"} else
        return("N")
  
  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  
  return(x)
}


get_weak1 <- function(x, type){
  temp <- ""
  if (type){
    target <- str_sub(x , 3, 3)
  } else
    target <- str_sub(x , -3, -3)
  
  if (target == "C") {temp <- "A"} else
    if (target == "A") {temp <- "C"} else
      if (target == "G") {temp <- "T"} else
        if (target == "T") {temp <- "G"}
  if(type){
    substring(x, 3, 3) <- temp
  }else
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
  return(x)
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

## Clean out the nested table for the algorithum
get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}


# Get all endpoints and give parents and children
get_endpoints <- function(lst, current_name = "", parent_names = character()) {
  endpoints <- list()
  
  if (is.list(lst)) {
    if (length(lst) > 0) {
      for (i in seq_along(lst)) {
        nested_name <- names(lst)[i]
        nested_value <- lst[[i]]
        
        if (is.list(nested_value)) {
          nested_endpoints <- get_endpoints(
            nested_value,
            paste(current_name, nested_name, sep = "/"),
            c(parent_names, current_name)
          )
          endpoints <- c(endpoints, nested_endpoints)
        } else {
          endpoint <- list(
            endpoint = nested_name,
            parents = c(parent_names, current_name)
          )
          endpoints <- c(endpoints, list(endpoint))
        }
      }
    } else {
      endpoint <- list(
        endpoint = current_name,
        parents = parent_names
      )
      endpoints <- c(endpoints, list(endpoint))
    }
  } else {
    endpoint <- list(
      endpoint = current_name,
      parents = parent_names
    )
    endpoints <- c(endpoints, list(endpoint))
  }
  return(endpoints)
}


# A function tha cleans
clean_endpoints <- function(endpoints){
  for (i in 1:length(endpoints)){
    endpoints[[i]]$parents <- endpoints[[i]]$parents[-1]
    for (j in 1:length(endpoints[[i]]$parents)){
      
      split_string <- strsplit(endpoints[[i]]$parents[j], "/")[[1]]
      desired_item <- split_string[length(split_string)]
      endpoints[[i]]$parents[j] <- desired_item
    } 
  }
  return(endpoints)
}


# find the bad nodes
compute_bad_nodes <- function(endpoints, threshold){
  blacklist <- list()
  
  for (i in 1:length(endpoints)){
    result = 0
    for (j in 1:length(endpoints[[i]]$parents)){
        result = result + (calculate_dimer(endpoints[[i]]$endpoint, 
                                           endpoints[[i]]$parents[j])$temp > threshold)
        # print(calculate_dimer(endpoints[[i]]$endpoint, endpoints[[i]]$parents[j])$temp)
    }
    blacklist <- c(blacklist, result)
  }
  
  bad_nodes = endpoints[blacklist == 1]
  return(bad_nodes)
}


### Remove unqualafied nods from the OG df
remove_list <- function(lst, path) {
  if (length(path) == 1) {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- NULL
    }
  } else {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- remove_list(lst[[path[[1]]]], path[-1])
      if (is.list(lst[[path[[1]]]]) && length(lst[[path[[1]]]]) == 0 && !any(names(lst[[path[[1]]]]))) {
        lst[[path[[1]]]] <- NULL
      }
    }
  }
  lst
}

## remove the trunk
remove_empty_lists <- function(lst) {
  if (is.list(lst)) {
    lst <- lapply(lst, remove_empty_lists)
    lst <- lst[lengths(lst) > 0]
  }
  lst
}


### Remove based on index
Iterate_remove <- function(level3,bad_nodes){
  for (i in 1:length(bad_nodes)){
    level3 = remove_list(level3, c(bad_nodes[[i]]$parents, bad_nodes[[i]]$endpoint))
  }
  return(level3)
}

# Gather incoming list
incoming_list <- function(arranged_list){
  level4 <- list()
  for (item in arranged_list) {
    # Create a sublist with the name as the item
    sublist <- list(1)
    names(sublist) <- item
    level4[item] <- sublist
  }
  return(level4)
}

# replacing the end nodes so we have a new list at the end
replace_end_nodes <- function(lst, replace_lst) {
  if (is.list(lst)) {
    if (length(lst) == 0) {
      return(replace_lst)
    } else {
      return(lapply(lst, function(x) replace_end_nodes(x, replace_lst)))
    }
  } else {
    return(replace_lst)
  }
}


## get the primers that is around the SNP and apply mismatches rules
extract_substrings <- function(string, center, start_distance , end_distance) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()
  
  for (item in string) {
    
    # right flanking
    for (distance in start_distance:end_distance) {
      sub <- substr(item, center+1, center+1 + distance)
      substrings_right <- c(substrings_right,
                            toupper(reverseComplement(get_strong1(sub,1))),
                            toupper(reverseComplement(get_strong2(sub,1))),
                            toupper(reverseComplement(get_medium1(sub,1))),
                            toupper(reverseComplement(get_weak1(sub,1))))
    }
    
    # Left flanking
    for (distance in start_distance:end_distance) {
      # print(substr(string, center -distance, center))
      sub <- substr(item, center - distance, center +1)
      substrings_left <- c(substrings_left, 
                           get_strong1(sub,0),
                           get_strong2(sub,0),
                           get_medium1(sub,0),
                           get_weak1(sub,0))
    }
    
    start_distance = start_distance + 1
    end_distance = end_distance + 1
    
  }
  
  # Return the extracted substrings
  return(list(left = substrings_left[!substrings_left %in% "N"], 
              right = substrings_right[!substrings_right %in% "N"]))
}

## Get the primers that are far away from the SNP
extract_substrings_far <- function(string, 
                                   center, 
                                   start_distance , 
                                   end_distance, 
                                   far, 
                                   shift) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()
  
  for (i in seq(far, far + shift)){
    # to Right
    for (distance in start_distance:end_distance) {
      # print(substr(string, far + center, far + center + distance))
      sub = substr(string, i + center, i + center + distance)
      substrings_right <- c(substrings_right, toupper(reverseComplement(sub)))
    }
    
    # Left flanking
    for (distance in start_distance:end_distance) {
      # print(substr(string, center - distance - far, center - far))
      substrings_left <- c(substrings_left, 
                           substr(string, 
                                  center - distance - i, 
                                  center - i))}
  }
  
  # Return the extracted substrings
  return(list(right = substrings_right, left = substrings_left))
}


## The one that produce all the primers
all_text_warngling <- function(snp_wrangled, 
                               start_distance, 
                               end_distance, 
                               center, 
                               far, 
                               shift){
  
  ## extrac the candidate from the left side (upstream) of the SNP
  ## extrac the candidate from the right side (downstream) of the SNP
  ## We are only getting the primers that are closed to the SNP for now
  grouped_sequences <- snp_wrangled %>%
    group_by(snpID) %>%
    summarize(sequence_list = list(sequence)) %>% 
    mutate(substrings = map(sequence_list, ~extract_substrings(.x, 
                                                               center, 
                                                               start_distance,
                                                               end_distance))) %>% 
    unnest(substrings)
  
  
  grouped_sequences_far <- snp_wrangled %>% 
    group_by(snpID) %>%
    slice(1:1)%>%
    ungroup() %>% 
    mutate(substrings = map(sequence,
                            ~extract_substrings_far(.x,
                                                    center, 
                                                    start_distance, 
                                                    end_distance,
                                                    far, 
                                                    shift))) %>% unnest(substrings)
  
  grouped_sequences$faraway <- grouped_sequences_far$substrings
  grouped_sequences <-  grouped_sequences[, -2]
  
  
  grouped_sequences$direction <- duplicated(grouped_sequences[[1]]) 
  
  grouped_sequences <- grouped_sequences %>%
    mutate(direction = ifelse(direction == FALSE, "LEFT", "RIGHT"))
  
  return(grouped_sequences)
}


# Apply all the filters before multiplexing
stage1_filter <- function(df, 
                          desired_tm, 
                          diff, 
                          Homodimer,
                          hairpin){
  df
  for (i in 1:length(df[[2]])){
    df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_homodimer)[2,]) < Homodimer]
    df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_hairpin)[2,]) < hairpin]
    df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) < desired_tm + diff]
    df[[2]][[i]] <- df[[2]][[i]][unlist(sapply(df[[2]][[i]], calculate_tm)) > desired_tm - diff]
    if (length(df[[2]][[i]]) == 0){
      length(df[[3]][[i]]) <- 0
    }
  }
  df
  
  for (i in 1:length(df[[3]])){
    if (length(df[[3]][[i]]) != 0){
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) > desired_tm - diff]
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_hairpin)[2,]) < hairpin]
      }
  }
  df
  
  for (i in 1:length(df[[3]])){
    if (length(df[[3]][[i]]) != 0){
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_homodimer)[2,]) < Homodimer]
      df[[3]][[i]] <- df[[3]][[i]][unlist(sapply(df[[3]][[i]], calculate_tm)) < desired_tm + diff]
    }
  }
  df
  for (i in length(df[[1]]):1){
    if (length(df[[2]][[i]]) == 0){
      df <- df[-i, ]
    }
  }
  df
  
  return(df)
}


### select the top n primers for multiplexing
extract_top_n <- function(nested_list, n) {
  modified_list <- lapply(nested_list, function(inner_list) {
    if (length(inner_list) >= n) {
      inner_list[1:n]
    } else {
      inner_list
    }
  })
  return(modified_list)
}

# This handle what part of the tree we want to show
get_display_tree <- function(level3, keep){
  endpoints <- get_endpoints(level3)
  
  # Endpoints come back a little messy
  endpoints <- clean_endpoints(endpoints)
  
  
  display_tree <- list()
  for (i in 1:keep){
    display_tree <- c(display_tree, list(unlist(endpoints[[i]])))
  }
  
  display_tree <- data.frame(display_tree)
  
  first_row <- display_tree[1, ]
  display_tree <- display_tree[-1, ]
  
  display_tree <- rbind(display_tree, first_row)
  
  colnames(display_tree) <- paste0("Option ", seq(1, keep))
  
  
  return(display_tree)
}



soulofmultiplex <- function(df, Heterodimer_tm){
  list_3 <- list()
  
  for (i in 1:length(df[[1]])){
    list_3 <- c(list_3, 
                list(unlist(df[[4]][[i]])), 
                list(unlist(df[[5]][[i]])))
  }
  
  # Arrange the list from small to big
  arranged_list <- list_3
  
  # Prepare the initial list for multiplexing
  level2 <- list()
  level3 <- list()
  level4 <- list()
  
  level2 <- incoming_list(arranged_list[[1]])
  
  level3 <- replace_end_nodes(incoming_list(arranged_list[[1]]),
                              incoming_list(arranged_list[[2]])
  )
  
  if (length(arranged_list) != 2) {
  
  level3 <- replace_end_nodes(level3,
                              incoming_list(arranged_list[[3]])
  )
  
  str(level3)
  # arranged_list
  # Running
  print(length(arranged_list))
  for (i in 4:length(arranged_list)){
    
    # Start a timer
    start_time <- Sys.time()
    
    # Get all the end points from the tree
    endpoints <- get_endpoints(level3)
    
    # Endpoints come back a little messy
    endpoints <- clean_endpoints(endpoints)
    print(paste("Start with ", length(endpoints)))
    
    # Evalauate all the ned points to its parents
    bad_nodes <- compute_bad_nodes(endpoints, Heterodimer_tm)
    print(paste("We are removing: ", length(bad_nodes)))
    
    
    # Remove bad nodes if there are any
    if (length(bad_nodes) != 0){
      level3 <- Iterate_remove(level3,bad_nodes)
      level3 <- remove_empty_lists(level3)
    }
    
    
    # If all nodes are bad, return NULL
    if (length(endpoints) == length(bad_nodes)){
      print("All nodes are removed during the process")
      return(NULL)
    }
    
    print(paste("After trimming: ", length(get_endpoints(level3))))
    
    # Stop adding list if we are at the last level
    if (1){
      level4 <- incoming_list(arranged_list[[i]])
      print(paste("New list: ", length(level4)))
      
      level3 <- replace_end_nodes(level3, level4)
      print(paste("level3 + level4: ", length(get_endpoints(level3))))
    }
    
    # Summarize results for this level
    print(paste("How far are we: ", i))
    print(paste("Time" , round(Sys.time() - start_time, 1)))
    print("--------------------------")
  }
  }
  
  level5 <- get_display_tree(level3, 3)
  repeated_list <- rep(df[[1]], each = 2)
  suffix <- c("_forward", "_reverse")
  modified_list <- paste0(repeated_list, suffix[rep(1:length(suffix), length.out = length(repeated_list))])
  rownames(level5) <- modified_list
  
  level5 <- rbind(level5, Tm = c(round(mean(sapply(level5[[1]], calculate_tm)), 2),
                            round(mean(sapply(level5[[2]], calculate_tm)), 2),
                            round(mean(sapply(level5[[3]], calculate_tm)),2)))

  
  return(level5)
  
}


get_tm_for_all_primers <- function(level5) {
  
  
  level5_with_tm_result <- data.frame(matrix(NA, nrow = nrow(level5), ncol = 0))
  
  # Apply the 'calculate_tm' function to each column of the dataframe
  for (i in seq_along(level5)) {
    # Calculate TM for the column and round the result
    tm_results <- round(calculate_tm(level5[[i]]), 2)
    
    # Combine the original column with the TM results
    combined <- data.frame(level5[[i]], tm_results)
    
    # Set the column names for the combined columns
    original_col_name <- names(level5)[i]
    names(combined) <- c(original_col_name, paste0(original_col_name, "_tm"))
    
    # Bind the new combined columns to the result dataframe
    level5_with_tm_result <- cbind(level5_with_tm_result, combined)
  }
  
  # Remove the first column if it contains only NA values from the placeholder creation
  level5_with_tm_result <- level5_with_tm_result[, colSums(is.na(level5_with_tm_result)) < nrow(level5_with_tm_result)]
  
}
