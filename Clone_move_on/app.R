



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
    textInput(inputId = "primer_list", label = "Enter Primers", value = "rs25 rs16944 rs1884 rs17287498"),
    numericInput(inputId = "primer_away", label = "Primer Distance (bp)", value = 476),
    sliderInput("primer_left_length", label = ("Forward (bp)"), min = 10,
                max = 40, value = c(15, 20)),
    sliderInput("primer_right_length", label = ("Reverse (bp)"), min = 10,
                max = 40, value = c(15, 20)),
    sliderInput("left_TM", "Left TM max", 1, 100, 70),
    sliderInput("right_TM", "Right TM max", 1, 100, 70),
    sliderInput("left_hair_TM", "Left hairpin TM max", 1, 100, 70),
    sliderInput("right_hair_TM", "Right hairpin TM max", 1, 100, 70),
    sliderInput("diff", "Max difference in TM", 1, 10, 5),
    sliderInput("Homodimer_left_dg", "Homodimer_left_dg", 1, 10, 5),
    sliderInput("Homodimer_right_dg", "Homodimer_right_dg", 1, 10, 5),
    sliderInput("Heterodimer_dg", "Heterodimer_dg", 1, 10, 5)
  ),
  
  dashboardBody(
    fluidRow(
      column(
        plotOutput("snippet1"),
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
    target <- complement(target)
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
    target <- complement(target)
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
    target <- complement(target)
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
    target <- complement(target)
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
  
  
  
  # 
  # primer <- "rs25 rs16944 rs1884 rs17287498"
  # primer_away <- 100
  # primer_min <- 20
  # primer_max <- 30
  # primer_left_min <- 18
  # primer_left_max <- 30
  # left_TM <- 70
  # right_TM <- 70
  # left_hair_TM <- 70
  # right_hair_TM <- 70
  # diff <- 5
  # Homodimer_left_dg <- 5
  # Homodimer_right_dg <- 5
  # Heterodimer_dg <- 5
  
  mart_api <- function(primer,
                       primer_away,
                       primer_min,
                       primer_max,
                       primer_left_min,
                       primer_left_max, 
                       left_TM, 
                       right_TM, 
                       left_hair_TM, 
                       right_hair_TM, 
                       diff,
                       Homodimer_left_dg, 
                       Homodimer_right_dg, 
                       Heterodimer_dg){
  
    print("Check point 1")
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
    
    snp_sequence_split <- as.data.frame(str_split(snp_sequence$`Variant sequences`, "%"))
    colnames(snp_sequence_split) <- snp_sequence$`Variant name`
    rownames(snp_sequence_split) <- c("upstream", "variants", "downstream")
    snps <- t(snp_sequence_split)
    snpsTibble <- as_tibble(snps, rownames = NA)
    vars_split <- str_split(snpsTibble$variants, "/")
    n <- 10
    vars_split_uniform <- as.data.frame(lapply(vars_split, `length<-`, n))
    # need to name things to keep that straight
    colnames(vars_split_uniform) <- snp_sequence$`Variant name`
    vars_split_transposed <- t(vars_split_uniform)
    variantsTibble <- as_tibble(vars_split_transposed, rownames = NA)
    varListFinal <- mutate(variantsTibble, snpsTibble)
    varListFinal['snpID'] <- snp_sequence$`Variant name`
    variantsTibbleFinal <- as_tibble(varListFinal, rownames = NA)
    variantsTibbleFinal2 <- pivot_longer(variantsTibbleFinal,
                                         cols = V1:V10,
                                         names_to = "variations",
                                         values_to = "observations")
    
    variantsTrimmed <- drop_na(variantsTibbleFinal2)
    # variantsTrimmed has everything we need to make the output for calling primer3
    # need to reorder then combine some columns
    # we want one string that is the upstream sequence, then the observation of the snp
    # then the downstream sequence all together in one string
    variantsTrimmed <- variantsTrimmed %>% relocate(snpID)
    variantsTrimmed <- variantsTrimmed %>% relocate(observations, .before = variations)
    variantsTrimmed <- variantsTrimmed %>% relocate(downstream, .before = variations)
    variantsTrimmed <- variantsTrimmed %>% unite("sequence", upstream:downstream, sep = "")
    # add columns for the substrings leading up to and including the variant site
    for (i in primer_left_min:primer_left_max) {
      colname <- paste0("left", i)
      variantsTrimmed <- variantsTrimmed %>%
        mutate(!!colname := str_sub(sequence, 501 - i, 501))
    }
    for (i in primer_min:primer_max) {
      colname <- paste0("right", 500 - primer_away -i)
      variantsTrimmed <- variantsTrimmed %>% mutate(!!colname := str_sub(sequence,
                                                                         500 - primer_away - i,
                                                                         500 - primer_away))
    }
    limit_left_start <- paste("left", primer_left_max, sep = "")
    limit_left_stop <- paste("left", primer_left_min, sep = "")
    limit_right_start <- paste("right", 500 - primer_away - primer_max, sep = "")
    limit_right_stop <- paste("right", 500 - primer_away - primer_min, sep = "")
    
    variantsTrimmed2 <- pivot_longer(variantsTrimmed,
                                     cols = limit_left_start:limit_left_stop,
                                     names_to = "Left_side",
                                     values_to = "leftPrimers")
    variantsTrimmed2 <- pivot_longer(variantsTrimmed2,
                                     cols = limit_right_start:limit_right_stop,
                                     names_to = "Right_side",
                                     values_to = "rightPrimers")
    variantsTrimmed2 <- variantsTrimmed2[c(1,4,6,5,7)]
    
    mismatch_list <- variantsTrimmed2 %>%
      mutate(strong_mismatch_1 = map(leftPrimers, get_strong1),
             strong_mismatch_2 = map(leftPrimers, get_strong2),
             Medium_mismatch = map(leftPrimers, get_medium1),
             Weak_mismatch = map(leftPrimers, get_weak1)) %>%
      pivot_longer(
        cols = c(strong_mismatch_1,
                 strong_mismatch_2,
                 Medium_mismatch,
                 Weak_mismatch),
        names_to = "Mismatch",
        values_to = "primer",
        values_drop_na = TRUE) %>%
      mutate(Identidy = paste(snpID, Left_side, Right_side, Mismatch,sep = " ")) %>%
      as.data.frame() %>%
      dplyr::select(c(8, 7, 5)) %>%
      mutate(rightPrimers = toupper(reverseComplement(rightPrimers)))
    
    
    
    
    print("Check point 3")
    df <- get_data(mismatch_list, left_TM, right_TM, left_hair_TM, right_hair_TM, diff,
                   Homodimer_left_dg, Homodimer_right_dg, Heterodimer_dg)
    
    return(df)
  }
  
  

  masterTable <- reactive(mart_api(input$primer_list,
                                   input$primer_away,
                                   input$primer_right_length[1],
                                   input$primer_right_length[2],
                                   input$primer_left_length[1],
                                   input$primer_left_length[2],
                                   input$left_TM, 
                                   input$right_TM, 
                                   input$left_hair_TM, 
                                   input$right_hair_TM, 
                                   input$diff,
                                   input$Homodimer_left_dg, 
                                   input$Homodimer_right_dg, 
                                   input$Heterodimer_dg))
  
  
  output$primer_table <- renderDataTable(masterTable()
  )
  
  
  
  output$snippet1 <- renderPlot({
    
    df = masterTable()
    
    options(repr.plot.width=100, repr.plot.height=8)
    
    p1 <- ggplot(df, aes(`TM_left (°C)`, `TM_right (°C)`)) +
      geom_hex() + 
      labs(title="TM") +
      theme_classic()
    
    p2 <- ggplot(df, aes(`Hairpin_left (°C)`, `Hairpin_right (°C)`)) +
      geom_hex() + 
      labs(title="Hairpin") +
      theme_classic()
    
    p3 <- ggplot(df, aes(`Homodimer_Left (kcal/mol)`, `Homodimer_Right (kcal/mol)`)) +
      geom_hex() + 
      labs(title="Homodimer") +
      theme_classic()
    
    p1 + p2 + p3

  })

  
}

# Run the application 
shinyApp(ui = ui, server = server)
