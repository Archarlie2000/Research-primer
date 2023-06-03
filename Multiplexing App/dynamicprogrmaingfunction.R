

top <- 10
threshold = 3
level = 2
final = matrix()


df <- df2 %>% 
  group_by(Identify) %>%
  slice(1:top)%>%
  ungroup()

nested_tables <- split(df, df$Identify)
le <- length(nested_tables)


get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}

multiplex(get_list(1,2), get_list(2,2))

list_1 = get_list(1,2)
list_2 = get_list(2,2)


multiplex = function(list_1, list_2){
  if (length(list_1) != 0 && length(list_2) != 0 && level <=  le){
    level = level + 1
    combination <- expand.grid(append(list_1, list_2))
    indices <- evaluation(combination)
  
    output1 <- append(list_1[indices], list_2[indices])
  
    output2 <- get_list(level,2)
  
    final = combination[indices, ]
    multiplex(output1, output2)
  
  return(indices)
  }
  else{
    return(NULL)
  }
}


evaluation <- function(combinations){
  num_rows <- nrow(combinations)
  num_cols <- ncol(combinations)
  permutation <- factorial(num_cols)/factorial(num_cols-2)/2
  result_matrix <- matrix(0, nrow = num_rows, ncol = permutation)

  
  
  # Iterate over each row of the matrix
  for (i in 1:num_rows) {
    row <- combinations[i, ]
    print(paste("i is", i))
    # Convert elements of the row to strings
    row <- lapply(row, as.character)
    clock = 0
    # Iterate over each item in the row
    for (j in 1:num_cols) {
      # Compare the item with every other item in the row
      # print(paste("j is", j))
      
      if (max(result_matrix[i, ]) <= threshold){
        for (k in 1:num_cols) {
          # print(paste("k is", k))
          if (max(result_matrix[i, ]) <= threshold){
            if (k > j) {
              clock = clock + 1 
              result_matrix[i, clock] <- calculate_dimer(row[[j]], row[[k]])$temp
            }
          }
        }
      }
    }
  }
  
  result_matrix <- as.data.frame(result_matrix)
  row_indices <- which(apply(result_matrix, 1, function(row) all(row < 5)))
  
  return(row_indices)
}
