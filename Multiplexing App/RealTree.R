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
level2
newList



# Add the sublist to the new list
newList[[item]] <- sublist


names(sublist) <- "GTAA"
