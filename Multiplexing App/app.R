



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
  dashboardHeader(title = "Multiplexing"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Analysis", tabName = "Analysis", icon = icon("th")),
      menuItem("Selection", tabName = "Selection", icon = icon("th")),
      textInput(inputId = "primer_list", label = "Enter SNP", value = "rs17025867, rs9939609, rs7903146, rs1121980, rs76141775"),
      numericInput(inputId = "shift", label = "Shift (bp)", value = 600),
      numericInput(inputId = "desired_tm", label = "desired_tm (°C)", value = 60),
      sliderInput("diff", "Max difference in TM", 1, 5, 2),
      sliderInput("Homodimer_tm", "Homodimer (°C)", -20, 5, 0),
      sliderInput("Heterodimer_tm", "Heterodimer (°C)", -20, 5, -5)
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
  
  # primer = "rs17025867, rs9939609, rs7903146, rs1121980, rs76141775"
  # primer = "rs17025867"
  # shift = 600
  # desired_tm = 60
  # diff = 2
  # Homodimer_tm = 0
  # Heterodimer_tm = -5

  ## The main function
  mart_api <- function(primer,
                       shift){
    
    center <- 800
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
    print("Primer generate")
    return(df)
  }
  
  get_filter <- function(df,
                         desired_tm,
                         diff,
                         Homodimer_tm) {
    
    print("R get filter activated")
    df <- stage1_filter(df, desired_tm, diff, Homodimer)
    df
    
    print("Filtered")
    
    
    # Count how many candidates there are for each candidates
    df <- df %>%
      mutate(substrings_count = lengths(substrings),
             faraway_count = lengths(faraway)) %>%
      relocate(snpID, substrings_count, faraway_count, everything())
    
    # Display the updated nested tibble
    return(df)
  }
  
  get_multiplex <- function(df,
                            Heterodimer_tm){
    
    
    ## I have left and right primers for each SNP. Sometimes, I only
    # Have left or right for some SNP. I need to arrange a list that 
    # contain all SNP primers. However, for each set, it can be either
    # from the left or from the right. I can make several tress.
    # But how do I arrange list from the df?
    #
    
    
    ###
    top <- 7
    
    print("Tree search")
    
    # Keep only certain amount of candidates
    df[[4]] <- extract_top_n(df[[4]], top)
    df[[5]] <- extract_top_n(df[[5]], top)
    
    # This is a bug. I have not yet to figure out how to grwow
    # multiple tree. So I just take whatever left and dispose 
    # The other flanking direactions
    df <- df %>% 
      group_by(snpID) %>% 
      slice(1)
    
    df
    # Prepare the general list of multiplexing
    list_3 <- list()
    
    for (i in 1:length(df[[1]])){
      list_3 <- c(list_3, 
                  list(unlist(df[[4]][[i]])), 
                  list(unlist(df[[5]][[i]])))
    }
    
    # Arrange the list from small to big
    arranged_list <- list_3
    
    # Prepare the initial list for multiplexing
    level2 <- list()
    level3 <- list()
    level4 <- list()
  
    level2 <- incoming_list(arranged_list[[1]])
    
    level3 <- replace_end_nodes(incoming_list(arranged_list[[1]]),
                                incoming_list(arranged_list[[2]])
    )
    
    level3 <- replace_end_nodes(level3,
                                incoming_list(arranged_list[[3]])
    )
    
    # str(level3)
    # arranged_list
    # Running
    print(length(arranged_list))
    for (i in 4:length(arranged_list)){
      
      # Start a timer
      start_time <- Sys.time()
      
      # Get all the end points from the tree
      endpoints <- get_endpoints(level3)
      
      # Endpoints come back a little messy
      endpoints <- clean_endpoints(endpoints)
      print(paste("Start with ", length(endpoints)))
      
      # Evalauate all the ned points to its parents
      bad_nodes <- compute_bad_nodes(endpoints, Heterodimer_tm)
      print(paste("We are removing: ", length(bad_nodes)))
      
      
      # Remove bad nodes if there are any
      if (length(bad_nodes) != 0){
      level3 <- Iterate_remove(level3,bad_nodes)
      level3 <- remove_empty_lists(level3)
      }
      
      
      # If all nodes are bad, return NULL
      if (length(endpoints) == length(bad_nodes)){
        print("All nodes are removed during the process")
        return(NULL)
      }
      
      print(paste("After trimming: ", length(get_endpoints(level3))))
      
      # Stop adding list if we are at the last level
      if (1){
        level4 <- incoming_list(arranged_list[[i]])
        print(paste("New list: ", length(level4)))
      
        level3 <- replace_end_nodes(level3, level4)
        print(paste("level3 + level4: ", length(get_endpoints(level3))))
      }
      
      # Summarize results for this level
      print(paste("How far are we: ", i))
      print(paste("Time" , round(Sys.time() - start_time, 1)))
      print("--------------------------")
    }
    
    # This handle what part of the tree we want to show
    level5 <- get_display_tree(level3, 3)
    
    return(level5)
  }
  
  
  # This one produce the true table we used
  masterTable <- reactive(get_filter(unfiltered(),
                                     input$desired_tm,
                                     input$diff,
                                     input$Homodimer_tm
                                     ))
  
  # This produced the raw table that has not beend filtered
  unfiltered <- reactive(mart_api(input$primer_list,
                                  input$shift))

  
  # This produce summary of primer generations
  output$primer_table <- renderDataTable(
    masterTable()[c(1,2,3)]
  )
  
  # output$primer_table <- renderDataTable(mtcars)
  
  
  # THis produces the result of multiplexing
  output$multiplex_table <- renderDataTable(get_multiplex(masterTable(),
                                                          input$Heterodimer_tm)
  )

  
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
