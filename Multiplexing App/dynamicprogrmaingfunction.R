

top <- 10
threshold = 20
level = 2
final = matrix()



df2 <- read.csv("metadata.csv")[,-1]


df2$direction <- sapply(strsplit(as.character(df2$Identity), " "), function(x) x[[2]])
df2$Identity <- sapply(strsplit(as.character(df2$Identity), " "), function(x) x[[1]])


df <- df2 %>% 
  group_by(Identity) %>%
  slice(1:top)%>%
  ungroup()

nested_tables <- split(df, df$Identity)
le <- length(nested_tables)


get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}

# multiplex(get_list(1,2), get_list(2,2))

list_1 = c(list(unique(get_list(1,2)[[1]])), 
           list(unique(get_list(1,3)[[1]])),
           list(unique(get_list(2,2)[[1]])),
           list(unique(get_list(2,3)[[1]])))

list_2 = c(list(unique(get_list(3,2)[[1]])), 
           list(unique(get_list(3,3)[[1]])),
           list(unique(get_list(4,2)[[1]])),
           list(unique(get_list(4,3)[[1]])))


multiplex = function(list_1, list_2){
  if (length(list_1) != 0 && length(list_2) != 0 && level <=  le){
    level = level + 1
    combination <- expand.grid(append(list_1, list_2))
    indices <- evaluation(combination)
  
    output1 <- append(list(list_1[[1]][indices[1][[1]]]),
                      list(list_2[[1]][indices[2][[1]]]))
  
    output2 <- get_list(level,2)
  
    final = combination[indices, ]
    multiplex(output1, output2)
  
  return(indices)
  }
  else{
    return(NULL)
  }
}


evaluation <- function(combination){
  num_rows <- nrow(combination)
  num_cols <- ncol(combination)
  permutation <- factorial(num_cols)/factorial(num_cols-2)/2
  result_matrix <- matrix(0, nrow = num_rows, ncol = permutation)

  
  
  # Iterate over each row of the matrix
  start_time <- Sys.time()
  for (i in 1:num_rows) {
    row <- combination[i, ]
    row <- lapply(row, as.character)
    clock = 0
    #print(paste("i is", i))
    
    for (j in 1:num_cols) {

      #print(paste("j is", j))
      if (max(result_matrix[i, ]) <= threshold){
        for (k in j:num_cols) {
          if (k>j){
            #print(paste("k is", k))
            clock = clock + 1
            #print(paste("Order ----------- ", clock))
            #print(row[[j]])
            #print(row[[k]])
            result_matrix[i, clock] <- calculate_dimer(row[[j]], row[[k]])$temp
          }
        }
      }
    }
  }
  end_time <- Sys.time()
  print(elapsed_time <- end_time - start_time)
  
  result_matrix <- as.data.frame(result_matrix)
  row_indices <- which(apply(result_matrix, 1, 
                             function(row) all(row != 0)))
  # indices_1 <- unique(ceiling(row_indices/len_1))
  # indices_2 <- unique(sapply(row_indices, 
  #                            function(x) if(x%%len_2 == 0){return(len_2)}
  #                            else
  #                              {return(x%%len_2)} ))
  # print(indices_1)
  # print(indices_2)
  return(row_indices)
}
