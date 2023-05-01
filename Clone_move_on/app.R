



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


# Deployment
library(shinydashboard)
library(reticulate)
library(shiny)
library(rsconnect)
library(shinydashboard)


options(repos = BiocManager::repositories())


# Define UI for application that draws a histogram
ui <- dashboardPage(
  
  dashboardHeader(title = "Primer Selection"),
  
  dashboardSidebar(
    textInput(inputId = "primer_list", label = "Enter SNP", value = "rs25 rs16944 rs1884 rs17287498"),
    numericInput(inputId = "primer_away", label = "Amplicant Length (bp)", value = 400),
    sliderInput("primer_left_length", label = "Forward (bp)", min = 10,
                max = 40, value = c(15, 20)),
    sliderInput("primer_right_length", label = "Reverse (bp)", min = 10,
                max = 40, value = c(15, 20)),
    sliderInput("left_TM", "Left TM max", 1, 100, 70),
    sliderInput("right_TM", "Right TM max", 1, 100, 70),
    sliderInput("left_hair_TM", "Left hairpin TM max", 1, 100, 70),
    sliderInput("right_hair_TM", "Right hairpin TM max", 1, 100, 70),
    sliderInput("diff", "Max difference in TM", 1, 10, 2),
    sliderInput("Homodimer_left_dg", "Homodimer_left_dg", 1, 10, 5),
    sliderInput("Homodimer_right_dg", "Homodimer_right_dg", 1, 10, 5),
    sliderInput("Heterodimer_dg", "Heterodimer_dg", 1, 10, 5)
  ),
  
  dashboardBody(
    fluidRow(
      column(
        plotlyOutput("snippet1"),
        DT::dataTableOutput(outputId = "primer_table"), width = 12
      )
    ),
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Cpnnect R to
  source_python("getdata.py")
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
  })
  
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
  
  left_flanking_get_strong1 <- function(x){
    temp <- ""
    target <- str_sub(x , 3, 3)

    if (target == "A") {temp <- "G"} else
      if (target == "G") {temp <- "A"} else
        if (target == "C") {temp <- "T"} else
          if (target == "T") {temp <- "C"}
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
  ## Mismatching on Ts
  get_strong2 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "T") {
      temp <- "T"
      substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
      return(x)}
    else
      return(NULL)
  }
  
  left_flanking_get_strong2 <- function(x){
    temp <- ""
    target <- str_sub(x , 3, 3)

    if (target == "T") {
      temp <- "T"
      substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
      return(x)}
    else
      return(NULL)
  }
  
  
  
  get_medium1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "A") {temp <- "A"} else
      if (target == "G") {temp <- "G"} else
        if (target == "C") {temp <- "C"} else
          return(NULL)
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
  
  left_flanking_get_medium1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "A") {temp <- "A"} else
      if (target == "G") {temp <- "G"} else
        if (target == "C") {temp <- "C"} else
          return(NULL)
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
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
  
  left_flanking_get_weak1 <- function(x){
    temp <- ""
    target <- str_sub(x , - 3, - 3)

    if (target == "C") {temp <- "A"} else
      if (target == "A") {temp <- "C"} else
        if (target == "G") {temp <- "T"} else
          if (target == "T") {temp <- "G"}
    substring(x, nchar(x) - 2, nchar(x) - 2) <- temp
    return(x)
  }
  
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
  
  
  
  primer <- "rs16944"
  primer_away <- 400
  primer_min <- 15
  primer_max <- 16
  primer_left_min <- 18
  primer_left_max <- 19
  left_TM <- 70
  right_TM <- 70
  left_hair_TM <- 70
  right_hair_TM <- 70
  diff <- 5
  Homodimer_left_dg <- 5
  Homodimer_right_dg <- 5
  Heterodimer_dg <- 5
  
  mart_api <- function(primer,
                       primer_away,
                       primer_min,
                       primer_max,
                       primer_left_min,
                       primer_left_max){
    
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
    
    # add columns for the substrings leading up to and including the variant site
    for (i in primer_left_min:primer_left_max) {
      colname <- paste0("left", i)
      variantsTrimmed <- variantsTrimmed %>%
        mutate(!!colname := str_sub(sequence, 501 - i, 501))
    }
    
    for (i in primer_min:primer_max) {
      colname <- paste0("right", 500 - primer_away - i)
      variantsTrimmed <- variantsTrimmed %>% mutate(!!colname := str_sub(sequence,
                                                                         500 - primer_away - i,
                                                                         500 - primer_away))
    }
    
    for (i in primer_left_min:primer_left_max) {
      colname <- paste0("(left_flanking)_left", i)
      variantsTrimmed <- variantsTrimmed %>%
        mutate(!!colname := str_sub(sequence, 501, 500 + i))
    }
    
    for (i in primer_min:primer_max) {
      colname <- paste0("(left_flanking)_right", 500 + primer_away + i)
      variantsTrimmed <- variantsTrimmed %>% mutate(!!colname := str_sub(sequence,
                                                                         500 + primer_away - i,
                                                                         500 + primer_away))
    }
    
    
    limit_left_start <- paste("left", primer_left_max, sep = "")
    limit_left_stop <- paste("left", primer_left_min, sep = "")
    limit_right_start <- paste("right", 500 - primer_away - primer_max, sep = "")
    limit_right_stop <- paste("right", 500 - primer_away - primer_min, sep = "")

    
    left_flanking_limit_left_start <- paste("(left_flanking)_left", primer_left_max, sep = "")
    left_flanking_limit_left_stop <- paste("(left_flanking)_left", primer_left_min, sep = "")
    left_flanking_limit_right_start <- paste("(left_flanking)_right", 500 + primer_away + primer_max, sep = "")
    left_flanking_limit_right_stop <- paste("(left_flanking)_right", 500 + primer_away + primer_min, sep = "")
    
    
    variantsTrimmed2 <- pivot_longer(variantsTrimmed,
                                     cols = limit_left_start:limit_left_stop,
                                     names_to = "Left_side",
                                     values_to = "leftPrimers")
    variantsTrimmed2 <- pivot_longer(variantsTrimmed2,
                                     cols = limit_right_start:limit_right_stop,
                                     names_to = "Right_side",
                                     values_to = "rightPrimers")
    
    
    variantsTrimmed2 <- pivot_longer(variantsTrimmed2,
                                     cols = left_flanking_limit_left_start:left_flanking_limit_left_stop,
                                     names_to = "left_flanking_Left_side",
                                     values_to = "left_flanking_leftPrimers")
    variantsTrimmed2 <- pivot_longer(variantsTrimmed2,
                                     cols = left_flanking_limit_right_start:left_flanking_limit_right_stop,
                                     names_to = "left_flanking_Right_side",
                                     values_to = "left_flanking_rightPrimers")
    
    
    
    #variantsTrimmed2 <- variantsTrimmed2[c(1,4,6,5,7)]
    
  
    
    vt_partition_1 <- cbind(variantsTrimmed2$snpID, 
                              variantsTrimmed2$Left_side,
                              variantsTrimmed2$leftPrimers,
                              variantsTrimmed2$Right_side,
                              variantsTrimmed2$rightPrimers,
                              "right"
                              )
    
    vt_partition_2 <- cbind(variantsTrimmed2$snpID,
                            variantsTrimmed2$left_flanking_Left_side,
                            variantsTrimmed2$left_flanking_leftPrimers,
                            variantsTrimmed2$left_flanking_Right_side,
                            variantsTrimmed2$left_flanking_rightPrimers,
                            "left")
    
    variantsTrimmed2 <- rbind(vt_partition_1, vt_partition_2) %>% data.frame()
    
    
    colnames(variantsTrimmed2) <- c("snp", 
                                   "forward_position",
                                   "forward_primer",
                                   "reversed_position",
                                   "reversed_primer",
                                   "flanking_direction")
    
    
    
    variantsTrimmed2$forward_position <-  gsub("[(left_flanking)_]", "",
                                               as.character(variantsTrimmed2$forward_position))
    variantsTrimmed2$reversed_position <-  gsub("[(left_flanking)_right]", "",
                                             as.character(variantsTrimmed2$reversed_position))
    
    mismatch_list <- variantsTrimmed2 %>%
      mutate(
        strong_mismatch_1 = case_when(
          flanking_direction == "right" ~ get_strong1(forward_primer[1]),
          flanking_direction == "left" ~ left_flanking_get_strong1(forward_primer[1]) ),
        strong_mismatch_2 = case_when(
          flanking_direction == "right" ~ get_strong2(forward_primer[1]),
          flanking_direction == "left" ~ left_flanking_get_strong2(forward_primer[1]) ),
        Medium_mismatch = case_when(
          flanking_direction == "right" ~ get_medium1(forward_primer[1]),
          flanking_direction == "left" ~ left_flanking_get_medium1(forward_primer[1]) ),
        Weak_mismatch = case_when(
          flanking_direction == "right" ~ get_weak1(forward_primer[1]),
          flanking_direction == "left" ~left_flanking_get_weak1(forward_primer[1]) ))
    
  
    
    # %>%
    #   pivot_longer(
    #     cols = c(strong_mismatch_1,
    #              strong_mismatch_2,
    #              Medium_mismatch,
    #              Weak_mismatch),
    #     names_to = "Mismatch",
    #     values_to = "primer",
    #     values_drop_na = TRUE)
    
    
    # %>%
    #   mutate(Identidy = paste(snpID, Left_side, Right_side, Mismatch,sep = " ")) %>%
    #   as.data.frame() %>%
    #   dplyr::select(c(8, 7, 5)) %>%
    #   mutate(rightPrimers = toupper(reverseComplement(rightPrimers)))
    
    
    
    
    print("Produced Mismatch list")
    df <- get_data(mismatch_list)
    print("Produced Munfilter list")
    
    return(df)
  }
  
  
  source_python("getfilter.py")
  masterTable <- reactive(get_filter(unfiltered(),
                                     input$left_TM, 
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
                                   input$primer_left_length[2]))
  
  
  output$primer_table <- renderDataTable(masterTable()
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
