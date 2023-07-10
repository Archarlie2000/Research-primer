extract_substrings <- function(string, center, start_distance , end_distance) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()
  
  for (item in string) {
    center = center + 1
    start_distance = start_distance + 1
    end_distance = end_distance + 1
    
    
    # right flanking
    for (distance in start_distance:end_distance) {
      sub <- substr(item, center, center + distance)
      substrings_right <- c(substrings_right,
                            toupper(reverseComplement(get_strong1(sub))),
                            toupper(reverseComplement(get_strong2(sub))),
                            toupper(reverseComplement(get_medium1(sub))),
                            toupper(reverseComplement(get_weak1(sub))))
    }
    
    # Left flanking
    for (distance in start_distance:end_distance) {
      # print(substr(string, center -distance, center))
      sub <- substr(item, center -distance, center)
      substrings_left <- c(substrings_left, 
                           get_strong1(sub),
                           get_strong2(sub),
                           get_medium1(sub),
                           get_weak1(sub))
    }
    
  }
  
  
  
  # Return the extracted substrings
  return(list(left = substrings_left[!substrings_left %in% "N"], 
              right = substrings_right[!substrings_right %in% "N"]))
}








extract_substrings_far <- function(string, 
                                   center, 
                                   start_distance , 
                                   end_distance, 
                                   far, 
                                   shift) {
  # Empty lists to store substrings
  substrings_left <- list()
  substrings_right <- list()

    for (i in seq(far, far + shift)){
      # to Right
      for (distance in start_distance:end_distance) {
        # print(substr(string, far + center, far + center + distance))
        sub = substr(string, i + center, i + center + distance)
        substrings_right <- c(substrings_right, toupper(reverseComplement(sub)))
      }
      
      # Left flanking
      for (distance in start_distance:end_distance) {
        # print(substr(string, center - distance - far, center - far))
        substrings_left <- c(substrings_left, 
                             substr(string, 
                                    center - distance - i, 
                                    center - i))}
    }
  
  # Return the extracted substrings
  return(list(left = substrings_left, right = substrings_right))
}



all_text_warngling <- function(snp_wrangled, 
                   start_distance, 
                   end_distance, 
                   center, 
                   far, 
                   shift){
  
  center = 500
  grouped_sequences <- snp_wrangled %>%
    group_by(snpID) %>%
    summarize(sequence_list = list(sequence)) %>% 
    mutate(substrings = map(sequence_list, ~extract_substrings(.x, 
                                                               center, 
                                                               start_distance,
                                                               end_distance))) %>% 
    unnest(substrings)


  grouped_sequences_far <- snp_wrangled %>% 
    group_by(snpID) %>%
    slice(1:1)%>%
    ungroup() %>% 
    mutate(substrings = map(sequence,
                            ~extract_substrings_far(.x,
                                                    center, 
                                                    start_distance, 
                                                    end_distance,
                                                    far, 
                                                    shift))) %>% unnest(substrings)
  grouped_sequences$faraway <- grouped_sequences_far$substrings
  grouped_sequences <-  grouped_sequences[, -2]
  return(grouped_sequences)
}





