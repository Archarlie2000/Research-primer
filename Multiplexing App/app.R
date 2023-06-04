



#################################################
# Library installation


# install.packages("DT")
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("stringi")

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
      
      textInput(inputId = "primer_list", label = "Enter SNP", value = "rs1121980, rs9939609, rs7903146, rs7903146"),
      numericInput(inputId = "primer_away", label = "Interval (bp)", value = 350),
      numericInput(inputId = "shift", label = "Shift (bp)", value = 0),
      sliderInput("primer_left_length", label = "Forward (bp)", min = 15,
                  max = 30, value = c(18, 25)),
      sliderInput("primer_right_length", label = "Reverse (bp)", min = 15,
                  max = 30, value = c(18, 25)),
      sliderInput("left_TM", label = "FW TM max (°C)", min = 30,
                  max = 75, value = c(55, 70)),
      sliderInput("right_TM", "RV TM max (°C)", 1, 100, 70),
      sliderInput("left_hair_TM", "FW hairpin max (°C)", 1, 100, 70),
      sliderInput("right_hair_TM", "RV hairpin max (°C)", 1, 100, 70),
      sliderInput("diff", "Max difference in TM", 1, 10, 2),
      sliderInput("Homodimer_left_dg", "FW Homodimer (°C)", 1, 70,40),
      sliderInput("Homodimer_right_dg", "RV Homodimer (°C)", 1, 70, 40),
      sliderInput("Heterodimer_dg", "Heterodimer (°C)", 1, 70, 40)
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
          
          column(3, verbatimTextOutput('x4')),
          # DT::dataTableOutput(outputId = "multiplex_table"),
                )
          )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


  
  ## This is random graph that is for prove of concept
  output$distPlot <- renderPlot(mtcars)
  

  
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
    
    
    #start_time <- Sys.time()  
    df2$`TM_left (°C)` <- sapply(df2$Forward, calculate_tm)
    df2 <- df2[df2$`TM_left (°C)` < left_TM_max, ]
    df2 <- df2[df2$`TM_left (°C)` > left_TM_min, ]
    
    df2$`TM_right (°C)` <- sapply(df2$Reversed, calculate_tm)
    df2 <- df2[df2$`TM_right (°C)` < right_TM, ]
    
    df2$`TM_Diff (°C)` <- abs(df2$`TM_left (°C)` - df2$`TM_right (°C)`)
    df2 <- df2[df2$`TM_Diff (°C)` < diff, ]
    
    df2$`Hairpin_left (°C)` <- sapply(df2$Forward, function(x) calculate_hairpin(x)$temp)
    df2 <- df2[df2$`Hairpin_left (°C)` < left_hair_TM, ]
    
    df2$`Hairpin_right (°C)` <- sapply(df2$Reversed, function(x) calculate_hairpin(x)$temp)
    df2 <- df2[df2$`Hairpin_right (°C)` < right_hair_TM, ]
    
    df2$`Homodimer_Left (°C)` <- sapply(df2$Forward, function(x) calculate_homodimer(x)$temp)
    df2 <- df2[df2$`Homodimer_Left (°C)` < Homodimer_left_dg, ]
    
    df2$`Homodimer_Right (°C)` <- sapply(df2$Reversed, function(x) calculate_homodimer(x)$temp)
    df2 <- df2[df2$`Homodimer_Right (°C)` < Homodimer_right_dg, ]

    for (i in 1:nrow(df2)){
      df2$`Heterodimer (°C)`[i] = calculate_dimer(df2$Forward[i], df2$Reversed[i])$temp
    }
    df2 <- df2[df2$`Heterodimer (°C)` < Heterodimer_dg, ]

    

    
    # df2 <- df2[df2$`TM_left (°C)` < left_TM_max, ]
    # df2 <- df2[df2$`TM_left (°C)` > left_TM_min, ]
    # df2 <- df2[df2$`TM_right (°C)` < right_TM, ]
    # df2 <- df2[df2$`TM_Diff (°C)` < diff, ]
    # df2 <- df2[df2$`Hairpin_left (°C)` < left_hair_TM, ]
    # df2 <- df2[df2$`Hairpin_right (°C)` < right_hair_TM, ]
    # df2 <- df2[df2$`Heterodimer (°C)` < Heterodimer_dg, ]
    # df2 <- df2[df2$`Homodimer_Left (°C)` < Homodimer_left_dg, ]
    # df2 <- df2[df2$`Homodimer_Right (°C)` < Homodimer_right_dg, ]



    colnames(df2) <- c("Identity",
                       "Forward (bp)",
                       "Reversed (bp)",
                       "TM_F (°C)",
                       "TM_R (°C)",
                       "TM diff (°C)",
                       "TM_F Hairpin (°C)",
                       "TM_R Hairpin (°C)",
                       "Homodimer_F (°C)",
                       "Homodimer_R (°C)",
                       "Heterodimer (°C)")

    
    #df2 <- df2[ c(1,2,3,4,) ]
    df2 <- df2 %>% 
      mutate_if(is.numeric, round, digits = 2) %>% 
      arrange('TM diff (°C)', 'TM_L Hairpin (°C)')
    print(nrow(df2))
    print("Give df2")
    
    # write.csv(df2, "metadata.csv")
    
    
    return(df2)
  }
  

  # These are the paramters used for trouble shotting

  # primer <- "rs1121980, rs9939609, rs7903146, rs7903146, rs4402960"
  # primer_away <- 450
  # primer_min <- 18
  # primer_max <- 25
  # primer_left_min <- 18
  # primer_left_max <- 25
  # left_TM <- 70
  # right_TM <- 70
  # left_hair_TM <- 35
  # right_hair_TM <- 35
  # diff <- 2
  # Homodimer_left_dg <- 30
  # Homodimer_right_dg <- 30
  # Heterodimer_dg <- 10
  # shift <- 150
  # left_TM_max = 68
  # left_TM_min = 60

  
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
    
    for (primer_away in primer_away:(primer_away + shift)){
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
                            variantsTrimmed_temp_2$left_flanking_Right_side,
                            variantsTrimmed_temp_2$left_flanking_rightPrimers,
                            variantsTrimmed_temp_2$left_flanking_Left_side,
                            variantsTrimmed_temp_2$left_flanking_leftPrimers,
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
      {variantsTrimmed2$reversed_primer[i] <-  toupper(reverseComplement(variantsTrimmed2$reversed_primer[i]))
        variantsTrimmed2$strong_mismatch_1[i] <-  get_strong1(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$strong_mismatch_2[i] <-  get_strong2(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$medium_mismatch[i] <-  get_medium1(variantsTrimmed2$forward_primer[i])
       variantsTrimmed2$weak_mismatch[i] <-  get_weak1(variantsTrimmed2$forward_primer[i])
  
      } else
      {
      variantsTrimmed2$forward_primer[i] <-  toupper(reverseComplement(variantsTrimmed2$forward_primer[i]))
      variantsTrimmed2$strong_mismatch_1[i] <-  left_flanking_get_strong1(variantsTrimmed2$forward_primer[i])
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
      mutate(identity = paste(snp, flanking_direction)) %>%
      as.data.frame() %>%
      dplyr::select(c(9, 8, 5))
     
    colnames(mismatch_list) = c("identity", "Forward", "Reversed")

    
    print("nrow(mismatch_list)")
    print(nrow(mismatch_list))
    
    
    mismatch_list_collected <- rbind(mismatch_list_collected, mismatch_list)
  
    }
    
    print("rows of mismatch collected = ")
    print(nrow(mismatch_list_collected))
    
    # df <- mismatch_list_collected
    
    return(mismatch_list_collected)
  }
  
  get_multiplex <- function(df,
                            left_TM,
                            diff,
                            hetero){
    
    
    
    s = input$primer_table_rows_selected
    
    df = df[s,]
    
    print("multiplex activated")
    outputframe <- data.frame(matrix(ncol = 3, nrow = nrow(df)*nrow(df)))
    colnames(outputframe) <- c("name", "forward", "reverse")


    # Match every forward with every reverse primer for cross checking
    k <- 0
    
    
    print("total number = ")
    print(nrow(df)*nrow(df))
    for (i in 1:nrow(df)){
      for (j in 1:nrow(df)){
        k <- j+(i-1)*nrow(df)
        outputframe[k,1] <- paste(df[i,1], df[j,1])
        outputframe[k,2] <- df[i,2]
        outputframe[k,3] <- df[j,3]
      }
    }

    df2 <- outputframe


      
      
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
  
  
  output$multiplex_table <- renderDataTable(get_multiplex(masterTable(),
                                                          input$left_TM[1],
                                                          input$diff,
                                                          input$Heterodimer_dg)
  )

  
  
  output$x4 = renderPrint({
    s = input$primer_table_rows_selected
    if (length(s)) {
      cat('These rows were selected:\n\n')
      cat(s, sep = ', ')
    }
  })
  
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
