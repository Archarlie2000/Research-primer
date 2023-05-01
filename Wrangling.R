primer <- "rs25 rs16944 rs1884 rs17287498"
snp_list <- strsplit(primer, " ")[[1]]
upStream <- c("500")
downStream <- c("500")




snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
                      filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
                      checkFilters = FALSE,
                      values = list(snp_list, upStream, downStream),
                      mart = snpmart,
                      bmHeader = TRUE)

list_seq <- function(snp) {
  first_position <- unlist(gregexpr('/', snp))[1]
  number_slash <- str_count(snp, "/")
  block <- str_sub(snp, first_position -1 , 
                   first_position - 2 + number_slash*2 + 1) 
  block <- gsub("/", "", block)
  block <- strsplit(block, "")
  
  for (i in block) {
  k <- paste(str_sub(snp, 1, first_position - 2),
              i,
              str_sub(snp, first_position - 2 + number_slash * 2 + 2, str_length(snp) ),
              sep = "")
  
  }
  k = k[!grepl("W", k)]
  return(k)
  
  }



snp_wrangled <- data.frame (mutation  = c(),
                         sequence = c()
)


snp_wrangled <- data.frame(matrix(ncol = 2, nrow = 0))


for (j in snp_sequence$`Variant name`){
for (i in list_seq(snp_sequence$`Variant sequences`[snp_sequence$`Variant name`==j])){
snp_wrangled[nrow(snp_wrangled) + 1,] <- c(j, i)
}
  
}

colnames(snp_wrangled) = c("name", "sequence")

