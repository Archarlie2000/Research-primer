



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
# install.packages("shinycssloaders")

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
library(shinycssloaders)
library( shinyWidgets)

source("functions.R")

options(repos = BiocManager::repositories())
options(scipen = 999)


ui <- dashboardPage(
  dashboardHeader(title = "Acorn Finder 12/24/2023"),
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      # actionButton("run_button", "Run Analysis", icon = icon("play")),
      # numericInput(inputId = "shift", label = "Max Length (bp)", value = 400),
      # numericInput(inputId = "desired_tm", label = "desired_tm (°C)", value = 60),
      # sliderInput("diff", "Max difference in TM", 1, 10, 5),
      numericInput(inputId = "Heterodimer_tm", label = "Heterodimer (°C)", value = 50),
      numericInput(inputId = "Homodimer", label = "Homodimer (°C)", value = 30),
      numericInput(inputId = "top", label = "Top", value = 2),
      numericInput(inputId = "hairpin", label = "Max Hairpin (°C)", value = 45),
      div(style = "display: none", downloadButton("downloadData", "Download"))
      
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "dashboard",
              # Introduction and Instructions at the top
              HTML('<h1 style="padding-left: 20px; margin: 5px; font-size: 80px;">Acorn Finder</h1>'),
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">Acorn Finder is a comprehensive tool designed for rapid Allel-specific Multiplexing Primer Generation</h2>'),
              HTML('<hr style="border-top: 1px solid #ccc; margin-top: 10px; margin-bottom: 10px;">'),
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">1. Enter SNP IDs</h2>'),
              div(style = "color: blue; gradient; padding-left: 30px;", textInput(inputId = "primer_list", label = "", value = "rs53576, rs1815739, rs7412, rs429358, rs6152")),
              # Existing Dashboard content
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">2. Specify Desired Tm (°C)</h2>'),
              div(style = "padding-left: 30px;", sliderInput(inputId = "desired_tm", label = "", value = 60, min = 40, max = 80)),
              
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">3. Specify Max Length (bp)</h2>'),
              div(style = "padding-left: 30px;", sliderInput(inputId = "shift", label = "", value = 100, min = 100, max = 500)),
              
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">4. Specify Max difference in TM (°C)</h2>'),
              div(style = "padding-left: 30px;", sliderInput(inputId = "diff", label = "", value = 5, min = 0, max = 10)),
              
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">5. Adjust other filters as needed</h2>'),
              div(style = "padding-left: 30px;", actionButton("run_button", "Run Analysis", icon = icon("play"))),
              
              # verbatimTextOutput("consoleOutput"),
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">6. This table shows how many primers we are using from each SNP to perform multiplexing</h2>'),
              div(style = "padding-left: 30px;", 
                  column(
                    withSpinner(DT::dataTableOutput(outputId = "primer_table")), 
                    width = 12
                  )
              ),
              HTML('<h2 style="padding-left: 20px; margin: 5px; font-size: 20px;">7. Final result</h2>'),
              div(style = "padding-left: 30px;", 
                  withSpinner(DT::dataTableOutput(outputId = "multiplex_table"))
              )
      )
      # Removed the About tab
    )
  )
)




# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # These are the paramters used for trouble shooting
# 
# primer = "rs53576, rs1815739, rs7412, rs429358, rs6152"
# shift = 100
# desired_tm = 64
# diff = 3
# Heterodimer_tm = 50
# Homodimer <- 45
# top <- 2
  consoleText <- reactiveVal("")
  
  appendConsole <- function(message) {
    consoleText(paste0(consoleText(), "\n", message))
  }

  observe({
    input$someActionButton
    isolate({
      current_time <- Sys.time()
      message <- paste("Button clicked at:", current_time)
      appendConsole(message)
    })
  })
  
  # Render the console text
  output$consoleOutput <- renderText({
    consoleText()
  })
  

  mart_api <- function(primer,
                       shift, appendConsole){

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
    
    appendConsole("Primer generated")
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
    # Techincal debt
    df <- df[!duplicated(df$snpID), ]
    
    
    
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
  unfiltered <- eventReactive(input$run_button, {
    mart_api(input$primer_list, input$shift, appendConsole)
  })
  
  
  
  
  
  # This produce summary of primer generations
  output$primer_table <- renderDataTable(
    masterTable()[c(1,2,3)]
  )

  
  # output$primer_table <- renderDataTable(mtcars)
  
  
  # This produces the result of multiplexing
  output$multiplex_table <- renderDataTable({
    req(masterTable())  # Ensure that masterTable is not NULL
    get_multiplex(masterTable(),
                  input$Heterodimer_tm,
                  input$top)
  })
  
  
  
  # Download dataframe
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     paste("Save a df", ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(masterTable(), file, row.names = FALSE)
  #   }
  # )
}

# Run the application 
shinyApp(ui = ui, server = server)