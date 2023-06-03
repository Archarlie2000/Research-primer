
top <- 5


df <- df2 %>% 
  group_by(Identify) %>%
  slice(1:top)%>%
  ungroup()

nested_tables <- split(df, df$Identify)


Tm_lists <- list()

for (i in 1:length(nested_tables)){
  k <- str_flatten(nested_tables[[i]][[4]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  
  j <- str_flatten(nested_tables[[i]][[5]], collapse = " ")
  j <- as.list(strsplit(j, " "))
  Tm_lists <- append(Tm_lists, k)
  Tm_lists <- append(Tm_lists, j)
}


all_lists <- list()

for (i in 1:length(nested_tables)){
  k <- str_flatten(nested_tables[[i]][[2]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  
  j <- str_flatten(nested_tables[[i]][[3]], collapse = " ")
  j <- as.list(strsplit(j, " "))
  all_lists <- append(all_lists, k)
  all_lists <- append(all_lists, j)
}


TmCombination  <- expand.grid(Tm_lists)
TmCombination_df <- as.data.frame(lapply(TmCombination, as.numeric))

combinations <- expand.grid(all_lists)
combinations_df <- as.data.frame(combinations)





#######################################


num_rows <- nrow(combinations)
num_cols <- ncol(combinations)
permutation <- factorial(num_cols)/factorial(num_cols-2)/2
result_matrix <- matrix(0, nrow = num_rows, ncol = permutation)
threshold = 3


# Iterate over each row of the matrix
for (i in 1:num_rows) {
  row <- combinations[i, ]
  
  # Convert elements of the row to strings
  row <- lapply(row, as.character)
  clock = 0
  print(i)
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
row_indices <- which(apply(result_matrix, 1, function(row) all(row < 0)))
combinations2 <- combinations[row_indices, ]


########################################




