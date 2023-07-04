# get_endpoints <- function(lst, current_name = "", parent_names = character()) {
#   endpoints <- list()
#   
#   if (is.list(lst)) {
#     if (length(lst) > 0) {
#       for (i in seq_along(lst)) {
#         nested_name <- names(lst)[i]
#         nested_value <- lst[[i]]
#         
#         if (is.list(nested_value)) {
#           nested_endpoints <- get_endpoints(
#             nested_value,
#             paste(current_name, nested_name, sep = "/"),
#             c(parent_names, current_name)
#           )
#           endpoints <- c(endpoints, nested_endpoints)
#         } else {
#           endpoint <- list(
#             endpoint = nested_name,
#             parents = c(parent_names, current_name)
#           )
#           endpoints <- c(endpoints, list(endpoint))
#         }
#       }
#     } else {
#       endpoint <- list(
#         endpoint = current_name,
#         parents = parent_names
#       )
#       endpoints <- c(endpoints, list(endpoint))
#     }
#   } else {
#     endpoint <- list(
#       endpoint = current_name,
#       parents = parent_names
#     )
#     endpoints <- c(endpoints, list(endpoint))
#   }
#   
#   return(endpoints)
# }




endpoints <- get_endpoints(level3)
print(endpoints)

endpoints <- clean_endpoints(endpoints)
print(endpoints)

# Clean the formatting


# print(endpoints)

#####################################


# compute blacklist

# compute_bad_nodes <- function(endpoints){
#   blacklist <- list()
#   
#   for (i in 1:length(endpoints)){
#     result = 0
#     for (j in 1:length(endpoints[[i]]$parents)){
#       result = result + (calculate_dimer(endpoints[[i]]$endpoint, endpoints[[i]]$parents[j])$temp > threshold)
#     }
#     blacklist <- c(blacklist, result)
#   }
#   
#   bad_nodes = endpoints[blacklist == 1]
#   return(bad_nodes)
# }

bad_nodes <- compute_bad_nodes(endpoints)

# for (i in 1:length(endpoints)){
#   result = 0
#   for (j in 1:length(endpoints[[i]]$parents)){
#    result = result + (calculate_dimer(endpoints[[i]]$endpoint, endpoints[[i]]$parents[j])$temp > threshold)
#   }
#   blacklist <- c(blacklist, result)
# }
# 
# blacklist
# 
# bad_nodes <- endpoints[blacklist == 1]
print(bad_nodes)

# Print the filtered list
# print(endpoints)
# str(endpoints)
  
#############################################################
# Removing bad nodes

# remove_list <- function(lst, path) {
#   if (length(path) == 1) {
#     if (is.list(lst) && path[[1]] %in% names(lst)) {
#       lst[[path[[1]]]] <- NULL
#     }
#   } else {
#     if (is.list(lst) && path[[1]] %in% names(lst)) {
#       lst[[path[[1]]]] <- remove_list(lst[[path[[1]]]], path[-1])
#       if (is.list(lst[[path[[1]]]]) && length(lst[[path[[1]]]]) == 0 && !any(names(lst[[path[[1]]]]))) {
#         lst[[path[[1]]]] <- NULL
#       }
#     }
#   }
#   lst
# }
# 
# 
# remove_empty_lists <- function(lst) {
#   if (is.list(lst)) {
#     lst <- lapply(lst, remove_empty_lists)
#     lst <- lst[lengths(lst) > 0]
#   }
#   lst
# }

### Remove based on index

# Iterate_remove <- function(level3,bad_nodes){
#   for (i in 1:length(bad_nodes)){
#     level3 = remove_list(level3, c(bad_nodes[[i]]$parents, bad_nodes[[i]]$endpoint))
#   }
#   return(level3)
# }

level3 <- Iterate_remove(level3,bad_nodes)
level3 <- remove_empty_lists(level3)


str(level3)
# endpoints <- get_endpoints(level3)
# print(endpoints)


# str(level3)

###################################################


# 
# level4 <- list()
# for (item in get_list(2,3)[[1]]) {
#   # Create a sublist with the name as the item
#   sublist <- list(1)
#   names(sublist) <- item
#   level4[item] <- sublist
# }
# 
# # Gather incoming list
# incoming_list <- function(level3, arranged_list){
#   level4 <- list()
#   for (item in arranged_list) {
#     # Create a sublist with the name as the item
#     sublist <- list(1)
#     names(sublist) <- item
#     level4[item] <- sublist
#   }
#   return(level4)
# }

level4 <- incoming_list(level3, arranged_list[[3]])


# replace_end_nodes <- function(lst, replace_lst) {
#   if (is.list(lst)) {
#     if (length(lst) == 0) {
#       return(replace_lst)
#     } else {
#       return(lapply(lst, function(x) replace_end_nodes(x, replace_lst)))
#     }
#   } else {
#     return(replace_lst)
#   }
# }

# Define the list to be replaced


# Apply the function to the nested list
level3 <- replace_end_nodes(level3, level4)

# Print the updated list
print(level3)
str(level3)

endpoints <- get_endpoints(level3)
length(endpoints)





