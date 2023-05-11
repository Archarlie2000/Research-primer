



#################################################
# Library installation



# install.packages("reticulate")
# install.packages("DT")
# install.packages("shiny")
# install.packages("spgs")
# install.packages("rsconnect")
# install.packages("dplyr")
# install.packages("stringi")
# install.packages("tidyverse")
# install.packages("shinydashboard")
# install.packages("hexbin")
# install.packages("patchwork")
# install.packages("plotly")
# install.packages("TmCalculator")



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")



# Data processing
library(DT)
library(dplyr)
library(tidyverse)
library(stringi)

#graphing
library(ggplot2)
library(hexbin)
library(patchwork)
library(plotly)



# Bioinformatics
library(biomaRt)
library(spgs)
library(TmCalculator)


# Deployment
library(shinydashboard)
library(reticulate)
library(shiny)
library(rsconnect)
library(shinydashboard)


options(repos = BiocManager::repositories())


ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Analysis", tabName = "Analysis", icon = icon("th")),
      
      textInput(inputId = "primer_list", label = "Enter SNP", value = "rs1121980, rs9939609, rs7903146, rs7903146"),
      numericInput(inputId = "primer_away", label = "Amplicant Length (bp)", value = 350),
      numericInput(inputId = "shift", label = "Shift (bp)", value = 1),
      sliderInput("primer_left_length", label = "Forward (bp)", min = 15,
                  max = 30, value = c(18, 25)),
      sliderInput("primer_right_length", label = "Reverse (bp)", min = 15,
                  max = 30, value = c(18, 25)),
      sliderInput("left_TM", label = "Left TM max (°C)", min = 30,
                  max = 75, value = c(55, 70)),
      sliderInput("right_TM", "Right TM max", 1, 100, 70),
      sliderInput("left_hair_TM", "Left hairpin max (°C)", 1, 100, 70),
      sliderInput("right_hair_TM", "Right hairpin max (°C)", 1, 100, 70),
      sliderInput("diff", "Max difference in TM", 1, 10, 2),
      sliderInput("Homodimer_left_dg", "Homodimer left (°C)", 1, 70,40),
      sliderInput("Homodimer_right_dg", "Homodimer right (°C)", 1, 70, 40),
      sliderInput("Heterodimer_dg", "Heterodimer (°C)", 1, 70, 40)
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              column(
                DT::dataTableOutput(outputId = "primer_table"), width = 12,
              )
      ),
      
      # Second tab content
      tabItem(tabName = "Analysis",
              downloadButton("downloadData", "Download")
      ) 
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Cpnnect R to python file
  source_python("getdata.py")

  
  ## This is random graph that is for prove of concept
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
  })
  
  
  ## Get strong mismatch for the last three bp
  get_strong1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "A") {temp <- "G"} else
      if (target == "G") {temp <- "A"} else
        if (target == "C") {temp <- "T"} else
          if (target == "T") {temp <- "C"}
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
  ## Get strong mismatch for the first three bp
  left_flanking_get_strong1 <- function(x){
    temp <- ""
    target <- str_sub(x , 3, 3)

    if (target == "A") {temp <- "G"} else
      if (target == "G") {temp <- "A"} else
        if (target == "C") {temp <- "T"} else
          if (target == "T") {temp <- "C"}
    substring(x, 3, 3) <- temp
    return(x)
  }
  
  ## Get strong mismatch for the last three bp, there are two types
  get_strong2 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "T") {
      temp <- "T"
      substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
      return(x)}
    else
      return("N")
  }
  
  ## Get strong mismatch for the first three bp
  left_flanking_get_strong2 <- function(x){
    temp <- ""
    target <- str_sub(x , 3, 3)
    if (target == "T") {
      temp <- "T"
      substring(x, 3, 3) <- temp
      return(x)}
    else
      return("N")
  }
  
  
  ## Get medium mismatch for the last three bp
  get_medium1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "A") {temp <- "A"} else
      if (target == "G") {temp <- "G"} else
        if (target == "C") {temp <- "C"} else
          return("N")
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
    
  }
  
  ## Get weak mismatch for the first three bp
  left_flanking_get_medium1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "A") {temp <- "A"} else
      if (target == "G") {temp <- "G"} else
        if (target == "C") {temp <- "C"} else
          return("N")
    substring(x, 3, 3) <- temp
    return(x)
  }
  
  
  ## Get weak mismatch for the last three bp
  get_weak1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "C") {temp <- "A"} else
      if (target == "A") {temp <- "C"} else
        if (target == "G") {temp <- "T"} else
          if (target == "T") {temp <- "G"}
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
  
  ## Get weak mismatch for the first three bp
  left_flanking_get_weak1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "C") {temp <- "A"} else
      if (target == "A") {temp <- "C"} else
        if (target == "G") {temp <- "T"} else
          if (target == "T") {temp <- "G"}
    substring(x, 3, 3) <- temp
    return(x)
  }
  
  
  ## Get reversed sequnce of a string (There is a beter way to do this, but)
  reverse_chars <- function(string)
  {
    # split string by characters
    string_split = strsplit(string, split = "")
    # reverse order
    rev_order = nchar(string):1
    # reversed characters
    reversed_chars = string_split[[1]][rev_order]
    # collapse reversed characters
    paste(reversed_chars, collapse = "")
  }
  
  
  ## Warngling SNP list to individual rows
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
    k = gsub("%", "", k)
    k = k[!grepl("W", k)]
    return(k)
  }
  
  
  get_filter <- function(df, 
                         left_TM_min,
                         left_TM_max,
                         right_TM, 
                         left_hair_TM, 
                         right_hair_TM, 
                         diff, 
                         Homodimer_left_dg, 
                         Homodimer_right_dg, 
                         Heterodimer_dg) {
    
    print("R get filter activated")
    df2 <- df
    
    df2 <- df2[df2$`TM_left (°C)` < left_TM_max, ]
    df2 <- df2[df2$`TM_left (°C)` > left_TM_min, ]
    df2 <- df2[df2$`TM_right (°C)` < right_TM, ]
    df2 <- df2[df2$`TM_Diff (°C)` < diff, ]
    df2 <- df2[df2$`Hairpin_left (°C)` < left_hair_TM, ]
    df2 <- df2[df2$`Hairpin_right (°C)` < right_hair_TM, ]
    df2 <- df2[df2$`Heterodimer (kcal/mol)` < Heterodimer_dg, ]
    #df2 <- df2[df2$`Heterodimer (kcal/mol)` > 0, ]
    df2 <- df2[df2$`Homodimer_Left (kcal/mol)` < Homodimer_left_dg, ]
    #df2 <- df2[df2$`Homodimer_Left (kcal/mol)` > 0, ]
    df2 <- df2[df2$`Homodimer_Right (kcal/mol)` < Homodimer_right_dg, ]
    #df2 <- df2[df2$`Homodimer_Right (kcal/mol)` > 0, ]
    
    
    colnames(df2) <- c("Identidy",
                       "Forward (bp)",
                       "Reversed (bp)",
                       "TM_L (°C)",
                       "TM_R (°C)",
                       "TM diff (°C)",
                       "TM_L Hairpin (°C)",
                       "TM_R Hairpin (°C)",
                       "Heterodimer (°C)",
                       "Homodimer_L (°C)",
                       "Homodimer_R (°C)")
    
    #df2 <- df2[ c(1,2,3,4,) ]
    
    print("Give df2")
    
    
    
    return(df2)
  }
  

  # These are the paramters used for trouble shotting
  # primer <- "rs1121980, rs9939609"
  # primer_away <- 400
  # primer_min <- 15
  # primer_max <- 18
  # primer_left_min <- 15
  # primer_left_max <- 18
  # left_TM <- 70
  # right_TM <- 70
  # left_hair_TM <- 70
  # right_hair_TM <- 70
  # diff <- 5
  # Homodimer_left_dg <- 5
  # Homodimer_right_dg <- 5
  # Heterodimer_dg <- 5
  # shift <- 5
  
  
  ## The main function
  mart_api <- function(primer,
                       primer_away,
                       primer_min,
                       primer_max,
                       primer_left_min,
                       primer_left_max,
                       diff,
                       shift){
    
    ## Not sure why, but it works
    primer_away <- -primer_away
    
    
    # Accessing database
    print("Execute MART API")
    snp_list <- strsplit(primer, " ")[[1]]
    upStream <- c("500")
    downStream <- c("500")
    snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
    snp_sequence <- getBM(attributes = c('refsnp_id', 'snp'),
                          filters = c('snp_filter', 'upstream_flank', 'downstream_flank'),
                          checkFilters = FALSE,
                          values = list(snp_list, upStream, downStream),
                          mart = snpmart,
                          bmHeader = TRUE)
    

    
    ### Wrangling dataframe
    
    print("Data gathered")
    #Create a new data frame
    snp_wrangled <- data.frame(matrix(ncol = 2, nrow = 0))
    
    
    # Add each variation as a new string into each row
    for (j in snp_sequence$`Variant name`){
      for (i in list_seq(snp_sequence$`Variant sequences`[snp_sequence$`Variant name`==j])){
        
        snp_wrangled[nrow(snp_wrangled) + 1,] <- c(j, i)
      }
      
    }
    
    # Rename columns and data frame
    colnames(snp_wrangled) = c("snpID", "sequence")
    variantsTrimmed <- snp_wrangled
    variantsTrimmed_ghost <- snp_wrangled
    
    
    mismatch_list_collected <- data.frame(Identify = c(),
                                          Forward = c(),
                                          Reversed = c()
    ) 
    
    ### Insert mega for loop
    ###
    ###
    ###
    
    
    print("Start big loop")
    
    for (primer_away in primer_away:(primer_away+shift)){
      print("Amplicant distance")
      print(primer_away)
    
    # add columns for the substrings leading up to and including the variant site
    # produce right flanking left primer
    for (i in primer_left_min:primer_left_max) {
      colname <- paste0("left", i)
      variantsTrimmed <- variantsTrimmed %>%
        mutate(!!colname := str_sub(sequence, 501 - i, 501))
    }
    
    
    # produce right flanking right primer
    for (i in primer_min:primer_max) {
      colname <- paste0("right", 500 - primer_away - i)
      variantsTrimmed <- variantsTrimmed %>% mutate(!!colname := str_sub(sequence,
                                                                         500 - primer_away - i,
                                                                         500 - primer_away))
    }
    
    # produce left flanking left primer
    for (i in primer_left_min:primer_left_max) {
      colname <- paste0("(left_flanking)_left", i)
      variantsTrimmed_ghost <- variantsTrimmed_ghost %>%
        mutate(!!colname := str_sub(sequence, 501, 500 + i))
    }
    
    
    # produce left flanking right primer
    for (i in primer_min:primer_max) {
      colname <- paste0("(left_flanking)_right", 500 + primer_away + i)
      variantsTrimmed_ghost <- variantsTrimmed_ghost %>% mutate(!!colname := str_sub(sequence,
                                                                         500 + primer_away - i,
                                                                         500 + primer_away))
    }
    
    ## Define the range of flanking for pivoting (right flanking)
    limit_left_start <- paste("left", primer_left_max, sep = "")
    limit_left_stop <- paste("left", primer_left_min, sep = "")
    limit_right_start <- paste("right", 500 - primer_away - primer_max, sep = "")
    limit_right_stop <- paste("right", 500 - primer_away - primer_min, sep = "")

    ## Define the range of flanking for pivoting (left flanking)
    left_flanking_limit_left_start <- paste("(left_flanking)_left", primer_left_max, sep = "")
    left_flanking_limit_left_stop <- paste("(left_flanking)_left", primer_left_min, sep = "")
    left_flanking_limit_right_start <- paste("(left_flanking)_right", 500 + primer_away + primer_max, sep = "")
    left_flanking_limit_right_stop <- paste("(left_flanking)_right", 500 + primer_away + primer_min, sep = "")
    
    
    
    
    ## Pivot the column into a long list
    variantsTrimmed_temp_1 <- pivot_longer(variantsTrimmed,
                                     cols = limit_left_start:limit_left_stop,
                                     names_to = "Left_side",
                                     values_to = "leftPrimers") %>% 
                            pivot_longer(
                                     cols = limit_right_start:limit_right_stop,
                                     names_to = "Right_side",
                                     values_to = "rightPrimers")
    
    
    variantsTrimmed_temp_2 <- pivot_longer(variantsTrimmed_ghost,
                                     cols = left_flanking_limit_left_start:left_flanking_limit_left_stop,
                                     names_to = "left_flanking_Left_side",
                                     values_to = "left_flanking_leftPrimers") %>% 
                            pivot_longer(
                                     cols = left_flanking_limit_right_start:left_flanking_limit_right_stop,
                                     names_to = "left_flanking_Right_side",
                                     values_to = "left_flanking_rightPrimers")
    
    
    variantsTrimmed2 <- cbind(variantsTrimmed$snpID,
                              variantsTrimmed$sequence,
                              variantsTrimmed_temp_1, 
                              variantsTrimmed_temp_2)
    
    colnames(variantsTrimmed2)[1] <- "snpID"
    colnames(variantsTrimmed2)[2] <- "sequence"
    
    
    
    
    ## combine left and flanking into a longer list since 
    ## previous one is not split in the right way
    
    vt_partition_1 <- cbind(variantsTrimmed_temp_1$snpID, 
                            variantsTrimmed_temp_1$Left_side,
                            variantsTrimmed_temp_1$leftPrimers,
                            variantsTrimmed_temp_1$Right_side,
                            variantsTrimmed_temp_1$rightPrimers,
                              "right"
                              )
    
    vt_partition_2 <- cbind(variantsTrimmed_temp_2$snpID,
                            variantsTrimmed_temp_2$left_flanking_Left_side,
                            variantsTrimmed_temp_2$left_flanking_leftPrimers,
                            variantsTrimmed_temp_2$left_flanking_Right_side,
                            variantsTrimmed_temp_2$left_flanking_rightPrimers,
                            "left")
    
    variantsTrimmed2 <- rbind(vt_partition_1,vt_partition_2) %>% data.frame()
    
    
    
    
    colnames(variantsTrimmed2) <- c("snp", 
                                   "forward_position",
                                   "forward_primer",
                                   "reversed_position",
                                   "reversed_primer",
                                   "flanking_direction")
    
    
    ## Fix the syntax for naming
    variantsTrimmed2$forward_position <-  gsub("[(left_flanking)_]", "",
                                               as.character(variantsTrimmed2$forward_position))
    variantsTrimmed2$reversed_position <-  gsub("[(left_flanking)_right]", "",
                                             as.character(variantsTrimmed2$reversed_position))
    
    
    ### Get mismatches for left primers depend on the flanking direaction
    for (i in 1:nrow(variantsTrimmed2)){
      if (variantsTrimmed2$flanking_direction[i] == "right")
      {variantsTrimmed2$strong_mismatch_1[i] <-  get_strong1(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$strong_mismatch_2[i] <-  get_strong2(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$medium_mismatch[i] <-  get_medium1(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$weak_mismatch[i] <-  get_weak1(variantsTrimmed2$forward_primer[i])
      } else
      {variantsTrimmed2$strong_mismatch_1[i] <-  left_flanking_get_strong1(variantsTrimmed2$forward_primer[i])
      variantsTrimmed2$strong_mismatch_2[i] <-  left_flanking_get_strong2(variantsTrimmed2$forward_primer[i])
      variantsTrimmed2$medium_mismatch[i] <-  left_flanking_get_medium1(variantsTrimmed2$forward_primer[i])
      variantsTrimmed2$weak_mismatch[i] <-  left_flanking_get_weak1(variantsTrimmed2$forward_primer[i])}
      }
    
    
    
    ## Pivot all mismtaches into a long list
    ## Remove all Ns from primer list
    ## get reversed complementrayr for right primer
    ## deselect some columns
    mismatch_list <- variantsTrimmed2 %>%
      pivot_longer(
        cols = c(strong_mismatch_1,
                 strong_mismatch_2,
                 medium_mismatch,
                 weak_mismatch),
        names_to = "Mismatch",
        values_to = "primer",
        values_drop_na = TRUE) %>% 
      filter(primer != "N") %>% 
      mutate(Identidy = paste(snp, 
                              forward_position, 
                              reversed_position, 
                              Mismatch, 
                              flanking_direction,sep = ", ")) %>%
      as.data.frame() %>%
      mutate(reversed_primer = toupper(reverseComplement(reversed_primer))) %>% 
      dplyr::select(c(9, 8, 5))
     
    colnames(mismatch_list) = c("Identify", "Forward", "Reversed")

    
  
    ## Pre filtering with TM
    con <- list()
    for (i in 1:nrow(mismatch_list)){

      k = (Tm_NN(mismatch_list$Forward[i], Na =50)[[1]] - 
        Tm_NN(mismatch_list$Reversed[i], Na =50)[[1]]) <= diff
      
      con <- append(con, k )
    }
    
    con <- as.logical(con)
    mismatch_list$temp <- con
    
    # Drop the temp colum
    mismatch_list <- mismatch_list[mismatch_list$temp == 1, ]  %>%  
      dplyr::select (-temp)
    
    print("nrow(mismatch_list)")
    print(nrow(mismatch_list))
    
    
    mismatch_list_collected <- rbind(mismatch_list_collected, mismatch_list)
  
    }
    
    print("rows of mismatch collected = ")
    print(nrow(mismatch_list_collected))
    
    
    df <- get_data(mismatch_list_collected)
    print("Unfiltered list Produced ")
    
    return(df)
  }
  
  
  
  masterTable <- reactive(get_filter(unfiltered(),
                                     input$left_TM[1],
                                     input$left_TM[2],
                                     input$right_TM, 
                                     input$left_hair_TM, 
                                     input$right_hair_TM, 
                                     input$diff,
                                     input$Homodimer_left_dg, 
                                     input$Homodimer_right_dg, 
                                     input$Heterodimer_dg))
  
  
  unfiltered <- reactive(mart_api(input$primer_list,
                                   input$primer_away,
                                   input$primer_right_length[1],
                                   input$primer_right_length[2],
                                   input$primer_left_length[1],
                                   input$primer_left_length[2],
                                  input$diff,
                                  input$shift))
  
  
  output$primer_table <- renderDataTable(masterTable()
  )
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Save a df", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(masterTable(), file, row.names = FALSE)
    }
  )
  
  
  output$snippet1 <- renderPlotly({
    
    df = masterTable()
    
    options(repr.plot.width=100, repr.plot.height=8)
    
    print("get graphes")
    p1 <- ggplotly(ggplot(df, aes(`TM_left (°C)`, `TM_right (°C)`)) +
      geom_hex() + 
      labs(title="TM") +
      theme_classic())
    
    p2 <- ggplotly(ggplot(df, aes(`Hairpin_left (°C)`, `Hairpin_right (°C)`)) +
      geom_hex() + 
      labs(title="Hairpin") +
      theme_classic())
    
    p3 <- ggplotly(ggplot(df, aes(`Homodimer_Left (kcal/mol)`, `Homodimer_Right (kcal/mol)`)) +
      geom_hex() + 
      labs(title="Homodimer") +
      theme_classic())
    
    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    
    p1 %>% layout(autosize = F, width = 500, height = 400, margin = m)

  })

  
}

# Run the application 
shinyApp(ui = ui, server = server)
