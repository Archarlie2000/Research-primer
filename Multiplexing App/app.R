



#################################################
# Library installation


# install.packages("DT")
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("stringi")
# install.packages("stringr")
# install.packages("mosaic")
# install.packages("purrr")

# install.packages("ggplot2")
# install.packages("hexbin")
# install.packages("patchwork")
# install.packages("plotly")

# install.packages("devtools")
# devtools::install_github("jensenlab/primer3")
# install.packages("TmCalculator")
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
# install.packages("spgs")

# install.packages("shiny")
# install.packages("rsconnect")
# install.packages("shinydashboard")

# you mush downgrade your dbplyr package to avoid conflict with biomart
# devtools::install_version("dbplyr", version = "2.3.4")

# Probe
library(rprimer)
#library(Biostrings)

# Data processing
library(DT)
library(dplyr)
library(tidyverse)
library(stringi)
library(stringr)
library(mosaic)
library(purrr)


#graphing
library(ggplot2)
library(hexbin)
library(patchwork)
library(plotly)


# Bioinformatics
library(biomaRt)
library(spgs)
library(primer3)


# Deployment
library(shinydashboard)
library(shiny)


source("functions.R")

options(repos = BiocManager::repositories())

ui <- dashboardPage(
  dashboardHeader(title = "Acorn Finder"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Analysis", tabName = "Analysis", icon = icon("th")),
      menuItem("Selection", tabName = "Selection", icon = icon("th")),
      textInput(inputId = "primer_list", label = "Enter SNP", value = "rs53576, rs1815739, rs7412, rs429358, rs6152"),
      numericInput(inputId = "shift", label = "Shift (bp)", value = 100),
      numericInput(inputId = "desired_tm", label = "desired_tm (°C)", value = 60),
      sliderInput("diff", "Max difference in TM", 1, 10, 5),
      numericInput(inputId = "Heterodimer_tm", label = "Heterodimer (°C)", value = 50),
      numericInput(inputId = "Homodimer", label = "Homodimer (°C)", value = 30),
      numericInput(inputId = "top", label = "Top", value = 2),
      numericInput(inputId = "hairpin", label = "hairpin (°C)", value = 45)
    )
  ),
  dashboardBody(
    
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              column(
                DT::dataTableOutput(outputId = "primer_table"), 
                width = 12)),
      # Second tab content
      tabItem(tabName = "Analysis",
              downloadButton("downloadData", "Download"))
    ),
    tabItem(tabName = "Selection",
            
            DT::dataTableOutput(outputId = "multiplex_table"),
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # These are the paramters used for trouble shooting

  # primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
  # shift = 100
  # desired_tm = 60
  # diff = 5
  # Heterodimer_tm = 50
  # Homodimer <- 45
  # top <- 2

  
  ## The main function - generates our primers
  mart_api <- function(primer,
                       shift){

    # We will start exploring options 800 bp away from the SNP location upstream and downstream    
    center <- 800
    hairpin <- 45
    # from that distance of 800, we will search the range from 600 to 1,000. (800+200 and 800-200)
    far <- 200
    start_distance <- 15
    end_distance <- 30
    
    # Accessing database
    print("Execute MART API")
    snp_list <- strsplit(primer, " ")[[1]]
    upStream <- center
    downStream <- center
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
    
    
    ### I have a long long string. I want to get the left 18~25 charactors and 
    # between 300 ~ 800 units away, I want another 18 ~ 25
    df <- all_text_warngling(snp_wrangled, 
                             start_distance, 
                             end_distance, 
                             center, 
                             far, 
                             shift)
    
    
    
    df
    
    print("Primer generated")
    return(df)
  }
  
  get_filter <- function(df, # primer
                         desired_tm,
                         diff, # max diff in tm
                         Heterodimer_tm, # should this be heterodimer_tm?
                         Homodimer,
                         hairpin) {
    
    print("R get filter activated")
    # Applied filters before multiplexing
    df <- stage1_filter(df, desired_tm, diff, Homodimer, hairpin)
    print(df)
    
    print("Filtered")
    
    
    # Count how many candidates there are for each primer group
    df <- df %>%
      mutate(substrings_count = lengths(substrings),
             faraway_count = lengths(faraway)) %>%
      relocate(snpID, substrings_count, faraway_count, everything()) # Moves a block of columns
    
    # Display the updated nested tibble
    return(df)
  }
  
  get_multiplex <- function(df,
                            Heterodimer_tm,
                            top){
    
    print("Tree search")
    df
    # Keep only certain amount of candidates
    df[[4]] <- extract_top_n(df[[4]], top)
    df[[5]] <- extract_top_n(df[[5]], top)
    
    df <- df %>%
      group_by(snpID) %>%
      filter(substrings_count == max(substrings_count))
    
    print(df)
    
    
    level5 <- soulofmultiplex(df, Heterodimer_tm)
    print(level5)
  
    
    level5_with_tm_result <- get_tm_for_all_primers(level5)

    
    return(level5_with_tm_result)
  }
  
  
  # This one produces the true table we used
  masterTable <- reactive(get_filter(unfiltered(),
                                     input$desired_tm,
                                     input$diff,
                                     input$Heterodimer_tm,
                                     input$Homodimer,
                                     input$hairpin
  ))
  
  # This produced the raw table that has not been filtered
  unfiltered <- reactive(mart_api(input$primer_list,
                                  input$shift))
  
  
  # This produce summary of primer generations
  output$primer_table <- renderDataTable(
    masterTable()[c(1,2,3)]
  )

  
  # output$primer_table <- renderDataTable(mtcars)
  
  
  # This produces the result of multiplexing
  output$multiplex_table <- renderDataTable(get_multiplex(masterTable(),
                                                          input$Heterodimer_tm,
                                                          input$top)
  )
  
  # output$multiplex_table <- renderDataTable(mtcars)
  
  
  # Download dataframe
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Save a df", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(masterTable(), file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)