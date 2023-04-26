snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
                      filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
                      checkFilters = FALSE,
                      values = list(snp_list, upStream, downStream),
                      mart = snpmart,
                      bmHeader = TRUE)



snp <- "AAAACA/T/GGAAAA"

snp_trim <- list()


list_seq <- function(snp) {
  first_position <- unlist(gregexpr('/', snp))[1]
  number_slash <- str_count(snp, "/")
  block <- str_sub(snp, first_position -1 , 
                   first_position - 2 + number_slash*2 + 1) 
  block <- gsub("/", "", block)
  block <- strsplit(block, "")
  
  k <-  paste(str_sub(snp, 1, first_position - 2),
              i,
              str_sub(snp, first_position - 2 + number_slash * 2 + 2, str_length(snp) ),
              sep = ""
  )
  
  }

list_seq(snp)

my_list <- list_seq(snp)
my_list


my_list[! my_list %in% c('w')]

