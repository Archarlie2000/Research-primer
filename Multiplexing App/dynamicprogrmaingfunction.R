

top <- 3
threshold =  10
level = 2
final = list()



df2 <- read.csv("metadata.csv")[,-1]


df2$direction <- sapply(strsplit(as.character(df2$Identity), " "), function(x) x[[2]])
df2$Identity <- sapply(strsplit(as.character(df2$Identity), " "), function(x) x[[1]])


df <- df2 %>% 
  group_by(Identity) %>%
  slice(1:top)%>%
  ungroup()

nested_tables <- split(df, df$Identity)
levels <- length(nested_tables) *2


get_list <- function(i, j){
  k <- str_flatten(nested_tables[[i]][[j]], collapse = " ")
  k <- as.list(strsplit(k, " "))
  return(k)
}



list_3 = c(list(unique(get_list(1,2)[[1]])), 
           list(unique(get_list(1,3)[[1]])),
           list(unique(get_list(2,2)[[1]])),
           list(unique(get_list(2,3)[[1]])),
           list(unique(get_list(3,2)[[1]])), 
           list(unique(get_list(3,3)[[1]])),
           list(unique(get_list(4,2)[[1]])),
           list(unique(get_list(4,3)[[1]])))


arranged_list <- list_3[order(sapply(list_3, length), 
                                     decreasing = FALSE)]


print(arranged_list)


list_1 <- list(arranged_list[[1]])
list_2 <- list(arranged_list[[2]])

len_1 <- length(list_1[[1]])
len_2 <- length(list_2[[1]])

level = 3

####
list_1 <- output1
list_2 <- output2
####


multiplex(list_1, list_2)

multiplex = function(list_1, list_2){
  if (length(list_1) != 0 && length(list_2) != 0 && level <=  le){
    level = level + 1
    len_1 <- length(list_1[[1]])
    #print(len_1)
    len_2 <- length(list_2[[1]])
    #print(len_2)
    combination <- expand.grid(append(list_1, list_2))
    print(paste("Total rows", level, nrow(combination)))
    indices <- evaluation_new(combination, len_1, len_2)
    #print(indices)
    if (level == 2){
    output1 <- append(list(unlist(list_1)[unlist(indices[1])]),
                      list(unlist(list_2)[unlist(indices[2])]))
    }
    else
    {
      output1 <- append(list_1,
                        list(unlist(list_2)[unlist(indices[2])]))
    }
  
    #print(output1)
    output2 <- list(arranged_list[[level]])
    #print(output2)
    
    final <- output1
    
    multiplex(output1, output2)
  
  }
  else{
    return(output1)
  }
}


######################
# A for loop?


arranged_list <- list_3[order(sapply(list_3, length), 
                              decreasing = FALSE)]


print(arranged_list)


list_1 <- list(arranged_list[[1]])
list_2 <- list(arranged_list[[2]])

len_1 <- length(list_1[[1]])
len_2 <- length(list_2[[1]])
level <- 3

for (level in 3:levels+1){
  if (length(list_1) != 0 && length(list_2) != 0){
    len_1 <- length(list_1)
    print(paste("length 1:", len_1))
    len_2 <- length(list_2[[1]])
    print(paste("length 2:", len_2))
    combination <- expand.grid(append(list_1, list_2))
    print(paste("Total rows", level, "-------", nrow(combination)))
    indices <- evaluation_new(combination, len_1, len_2)
    # print(indices)
    if (level == 3){
      output1 <- append(list(unlist(list_1)[unlist(indices[1])]),
                        list(unlist(list_2)[unlist(indices[2])]))
    }
    else
    {
      output1 <- append(list_1,
                        list(unlist(list_2)[unlist(indices[2])]))
    }
    
    
    if (level == levels+1){
    
      final <- output1
    }
    else{
      # print("output1")
      # print(output1)
      output2 <- list(arranged_list[[level]])
      # print("output2")
      # print(output2)
      
      list_1 <- output1
      list_2 <- output2  
      
      
      final <- output1
    }
  }
}




