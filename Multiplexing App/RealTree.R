get_list(1,2)


newList <- list()

# Iterate over each item in the original list
for (item in get_list(1,2)[[1]]) {
  # Create a sublist with the name as the item
  sublist <- list(0)
  names(sublist) <- item
  newList[item] <- sublist
}

newList

level2 <- list()

for (item in get_list(1,3)[[1]]) {
  # Create a sublist with the name as the item
  sublist <- list(0)
  names(sublist) <- item
  level2[item] <- sublist

}
level2

level3 <- list()

for (i in names(newList)){
  level3[i] <- list(level2)
}


level3
# level2
# newList
# level3

level4 <- list()

for (item in get_list(2,2)[[1]]) {
  # Create a sublist with the name as the item
  sublist <- list(1)
  names(sublist) <- item
  level4[item] <- sublist
  }

level4
# level3
# names(level3)
# names(level3[[1]])

# level5 <- list()


 # level3[[1]][[1]] <- level4


for (j in names(level3)){
  for (k in names(level3[[j]])){
    level3[[j]][[k]] <- level4
  }
}
level3
str(level3)

# for (k in names(level3)){
#   for (i in names(level3[[k]])){
#     level5[k][i] <- list(level4)
#   }
# }

level5




# compare the name of the current list at the end point 
# to the parent node

names(level3[[1]])[[1]]

str(names(level3[[1]])[[1]])

class(names(level3[[1]])[[1]])

calculate_homodimer(names(level3[[1]])[[1]])$temp

names(level3)[[1]]

# Remove end point
level3[[1]][[2]][3] <- NULL


calculate_dimer(
  names(level3)[[1]],
  names(level3[[1]])[[1]]
  )$temp



# This function will compare everything come in to the last item

compare <- function(s1,s2,s3){
  
  print(paste("comparision 2: ", calculate_dimer(s1, s3)$temp
  ))
  print(paste("Comparision 1: ", calculate_dimer(s2, s3)$temp
  ))
  result = ( calculate_dimer(s1, s3)$temp < threshold) +
           ( calculate_dimer(s2, s3)$temp < threshold)
  
  return( (result == 0) )
}



# Evaluate parents
for (j in 1:length(names(level3))){
  for (k in 1:length(names(level3[[j]]))){
    for (m in 1:length(names(level3[[j]][[k]]))){
      
      
      print(paste("level 1: ", names(level3)[[j]], " ",j
                  ))
      print(paste("level 2: ", names(level3[[j]])[[k]], " ",k
                  ))
      print(paste("level 3: ", names(level3[[j]][[k]])[[m]], " ",m
                  ))
      
      s1 <- names(level3)[[j]]
      s2 <- names(level3[[j]])[[k]]
      s3 <- names(level3[[j]][[k]])[[m]]
      
      if (compare(s1,s2,s3)){
        print("not valid")
      }
      else{
        print("valid")
      }
      print("------------")
    }
  }
}






