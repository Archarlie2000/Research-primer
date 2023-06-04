
library(parallel)
library(doParallel)

#########################



numCores <- detectCores()
cl <- makeCluster(numCores)

registerDoParallel(cl)

clusterEvalQ(cl, {
  library(primer3)
  library(stringr)
})

###################
evaluation <- function(combination){
  num_rows <- nrow(combination)
  num_cols <- ncol(combination)
  permutation <- factorial(num_cols)/factorial(num_cols-2)/2
  result_matrix <- matrix(0, nrow = num_rows, ncol = permutation)
  
  
  # Iterate over each row of the matrix
  for (i in 1:num_rows) {
    row <- combination[i, ]
    row <- lapply(row, as.character)
    clock = 0
    print(paste("i is", i))
    
    for (j in 1:num_cols) {
      
      #print(paste("j is", j))
      if (1){
        for (k in j:num_cols) {
          if (k>j){
            clock = clock + 1
            result_matrix[i, clock] <- calculate_dimer(row[[j]], row[[k]])$temp
          }
        }
      }
    }
  }
  
  result_matrix <- as.data.frame(result_matrix)
  row_indices <- which(apply(result_matrix, 1, function(row) all(row < 5)))
  # indices_1 <- unique(ceiling(row_indices/len_1))
  # indices_2 <- unique(sapply(row_indices, 
  #                            function(x) if(x%%len_2 == 0){return(len_2)}
  #                            else
  #                              {return(x%%len_2)} ))
  # print(indices_1)
  # print(indices_2)
  return(result_matrix)
}


######################################################
#Partition

matrix_data <- combination

n <- 10

# Calculate the number of rows in each part
rows_per_part <- ceiling(nrow(matrix_data) / n)

# Divide the matrix into parts and store them in a list
divided_matrix <- split(matrix_data, rep(1:n, each = rows_per_part, length.out = nrow(matrix_data)))

# Iterate through the divided parts in the list
for (i in 1:n) {
  current_part <- divided_matrix[[i]]
  # Do something with the current_part
  print(current_part)
}



################################
{
# Apply the function in parallel

  start_time <- Sys.time()
  
  
  
  results <- foreach(i = divided_matrix) %dopar% {
    evaluation(i)
  }
  
  combined_matrix <- do.call(rbind, results)
  
  stopCluster(cl)
  #combined_results <- unlist(results)

  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)

}


