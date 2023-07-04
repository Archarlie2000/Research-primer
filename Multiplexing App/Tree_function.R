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
compute_bad_nodes <- function(endpoints){
  blacklist <- list()
  
  for (i in 1:length(endpoints)){
    result = 0
    for (j in 1:length(endpoints[[i]]$parents)){
      result = result + (calculate_dimer(endpoints[[i]]$endpoint, endpoints[[i]]$parents[j])$temp > threshold)
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
incoming_list <- function(level3, arranged_list){
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




