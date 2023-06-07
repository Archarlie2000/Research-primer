

top <- 10
threshold = 20
level = 2
final = matrix()
valid_option = 0



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


list_1 <- list(arranged_list[[1]])
list_2 <- list(arranged_list[[2]])

len_1 <- length(list_1[[1]])
len_2 <- length(list_2[[1]])


####
list_1 <- output1
list_2 <- output2
####

multiplex = function(list_1, list_2){
  if (length(list_1) != 0 && length(list_2) != 0 && level <=  le){
    level = level + 1
    len_1 <- length(list_1[[1]])
    print(len_1)
    len_2 <- length(list_2[[1]])
    print(len_2)
    combination <- expand.grid(append(list_1, list_2))
    indices <- evaluation_new(combination, len_1, len_2)
    print(indices)
    if (level == 2){
    output1 <- append(list(unlist(list_1)[unlist(indices[1])]),
                      list(unlist(list_2)[unlist(indices[2])]))
    }
    else
    {
      output1 <- append(list_1,
                        list(unlist(list_2)[unlist(indices[2])]))
    }
  
    print(output1)
    output2 <- list(arranged_list[[level]])
    print(output2)
    
    
    multiplex(output1, output2)
  
  return(indices)
  }
  else{
    return(NULL)
  }
}





