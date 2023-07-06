my_list <- list(
  "Level_A" = list(
    "Level_B" = list(
      "Level_C" = list(1,2,3),
      "Level_D" = 2
    ),
    "Level_E" = list(
      "Level_F" = 3,
      "Level_G" = 4,
      "Level_H" = 5
    )
  ),
  "Level_I" = list(
    "Level_J" = 6,
    "Level_K" = list(
      "Level_L" = 7,
      "Level_M" = 8,
      "Level_N" = list(
        "Level_O" = 9,
        "Level_P" = 10
      )
    )
  )
)

remove_item <- function(lst, path, item) {
  if (length(path) == 1) {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- Filter(function(x) x != item, lst[[path[[1]]]])
    }
  } else {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- remove_item(lst[[path[[1]]]], path[-1], item)
    }
  }
  lst
}


my_list <- remove_item(my_list, c("Level_A", "Level_B", "Level_C"), 2)

print(my_list)



############################################
remove_list <- function(lst, path) {
  if (length(path) == 1) {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- NULL
    }
  } else {
    if (is.list(lst) && path[[1]] %in% names(lst)) {
      lst[[path[[1]]]] <- remove_list(lst[[path[[1]]]], path[-1])
    }
  }
  lst
}

my_list <- list(
  "Level_A" = list(
    "level_B" = list(
      "level_C" = list(1, 2, 3),
      "level_D" = 2
    ),
    "level_E" = list(
      "level_F" = 3,
      "level_G" = 4,
      "level_H" = 5
    )
  ),
  "Level_I" = list(
    "level_J" = 6,
    "level_K" = list(
      "level_L" = 7,
      "level_M" = 8,
      "level_N" = list(
      )
    )
  )
)


my_list

my_list <- remove_list(my_list, c("Level_I", "level_K", "level_M"))

print(my_list)

my_list <- remove_list(level3, c("GTGCCTTCCTATGGTGCAGAGGATGA", 
                                 "GTGAGATTTCAGATCCACCTGCCTA", "GACGGGCCGTAACTGTTCCTCTCAAG"))

#########################
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

my_list <- list(
  "Level_A" = list(
    "level_B" = list(
      "level_C" = list(1, 2, 3),
      "level_D" = 2
    ),
    "level_E" = list(
      "level_F" = 3,
      "level_G" = 4,
      "level_H" = 5
    )
  ),
  "Level_I" = list(
    "level_J" = 6,
    "level_K" = list(
      "level_L" = 7,
      "level_M" = 8,
      "level_N" = list()
    )
  )
)

my_list <- remove_list(my_list, c("Level_A", "level_B", "level_C"))

print(my_list)


###################################
remove_empty_lists <- function(lst) {
  if (is.list(lst)) {
    lst <- lapply(lst, remove_empty_lists)
    lst <- lst[lengths(lst) > 0]
  }
  lst
}

my_list <- remove_empty_lists(my_list)


my_list

##################################################


append_list_to_end_nodes <- function(lst, append_lst) {
  if (is.list(lst)) {
    if (length(lst) == 0) {
      return(append_lst)
    } else {
      return(lapply(lst, function(x) append_list_to_end_nodes(x, append_lst)))
    }
  } else {
    return(list(lst, append_lst))
  }
}

# Define the list to be appended

# Apply the function to the nested list
updated_list <- append_list_to_end_nodes(level3, get_list(2,3))

# Print the updated list
print(updated_list)

str(updated_list)

######################################
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

# Define the list to be replaced
replace_lst <- level4

# Apply the function to the nested list
updated_list <- replace_end_nodes(level3, replace_lst)

# Print the updated list
print(updated_list)
str(updated_list)



level4 <- list()
for (item in get_list(2,3)[[1]]) {
  # Create a sublist with the name as the item
  sublist <- list(1)
  names(sublist) <- item
  level4[item] <- sublist
}


################################
