library(shiny)
library(shinythemes)
library(stringr)
library(DT)
library(dplyr)
library(ggforce)


source("code/least_squares_fit_model_shiny.R")

# Define the UI
ui <- fluidPage(
  #shinythemes::themeSelector(), 
  theme = shinytheme("united"),
  titlePanel(
    windowTitle = "Oskar's Tox Explorer!",
    div(
      style = "background-color: green; padding: 10px;",
      tags$h1("Ecotoxicity Data Explorer", style = "color: #333333;")
    )),
  
  tabsetPanel(
  
    tabPanel("SSD curves",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 selectInput("casNumberInput", "Select CAS Number Column:", choices = NULL, selectize = FALSE),
                 br(),
                 sliderInput("hcInput", "Select Response Level:", min = 5, max = 95, value = 20, step = 5, ticks = TRUE),
                 br(),
                 sliderInput("MC_n", "n Monte Carlo runs", min = 100000, max = 1000000, value = 100000, step = 100000, ticks = TRUE),
                 br(),
                 actionButton("runButton", "Run nls-model"),
                 br()
                 ),
               
               mainPanel(
                 width = 10,
                 br(), 
                 "Data summary of the nonlinear least squares model.",
                 dataTableOutput("outputNLSData", width = "auto"),
                 fluidRow(
                   column(width = 6,
                          br(), 
                          plotOutput("outputPlot", width = "100%", height = "600px")
                   ),
                   column(width = 6,
                          br(),
                          plotOutput("outputHist", width = "100%", height = "600px")
                   )
                 ),
                 br(),
                 "The raw data from which the model is based on.",
                 br(),
                 # Download dataset Button
                 "Download active dataset",
                 br(),
                 tags$style(type="text/css", "#downloadData {background-color: red; color: black;font-family: Courier New}"),
                 downloadButton("downloadData", "Download"),
                 br(),
                 #"A dataset showing the raw input data"
                 dataTableOutput("outputRAWData", width = "auto")
                            )
                   
             ) 
            ),
    
    tabPanel("HC20 uncertainty",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 class = "custom-sidebar",
                 textInput("search_2", "Search Value:", placeholder = "Enter search value"),
                 actionButton("filterButton_2", "Apply Filter"),
                 actionButton("clearButton_2", "Clear Filter"),
                 checkboxInput("OnlyValid", "Only view valid model data", value = FALSE)
               ),
               
               mainPanel(
                 width = 10,
                 div(style = "width:100%; overflow-x:auto",
                     DT::dataTableOutput("resultsTable_uncert")
                 )
               )
             ),
             
    ),
    
    tabPanel("HESTIA Ecotox database explorer",
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 class = "custom-sidebar",
                 #selectInput("column", "Select Column:", choices = c("CAS.Number", "PesticideAI_name", "Substance_type")),
                 #textInput("search", "Search Value:", placeholder = "Enter search value"),
                 #actionButton("filterButton_1", "Apply Filter"),
                 #actionButton("clearButton_1", "Clear Filter"),
                 #br(),
                 #br(),
                 br(),
                 checkboxGroupInput("Chem_categories", "Select Chemical Category:", 
                                    choices = c("Antibiotic", "Antiviral", "Other inorganic chemicals",
                                                "Other organic chemicals", "PPCP", "Pesticide",
                                                "Pharmaceutical", "Unknown"),
                                    selected = c("Antibiotic", "Antiviral", "Other inorganic chemicals",
                                                 "Other organic chemicals", "PPCP", "Pesticide",
                                                 "Pharmaceutical", "Unknown"))
               ),
               
               mainPanel(
                 width = 10,
                 br(),
                 br(),
                 br(),
                 div(style = "width:99%; height:100%; overflow-x:auto; overflow-y:auto;",
                     DT::dataTableOutput("resultsTable")
                 )
               )
             )
    )
  )
)


# Define the server logic
server <- function(input, output, session) {
  
  # Load the USETOX_adapted data set
  data <- reactive({
    read.csv("data/HESTIA_envirotox_cfs.csv") %>% 
      select(-definition)
  })
  # Load the NLS_output dataframe
  data_uncert <- reactive({
    read.csv("data/nls_output_df.csv")
  })
  
### Ecotox-database ### ------------------------------------------------
  # Filter the data based on user selections and display the results
  filteredData <- reactive({
    req(data())
    req(input$Chem_categories)
  #   
  #   column <- input$column
  #   searchValue <- input$search
  #   
  #   if (!is.null(column) && !is.null(searchValue)) {
  #     filteredData <- data()[grepl(searchValue, data()[[column]], ignore.case = TRUE), ]
  #   } else {
  #     filteredData <- data()
  #   }
    if (length(input$Chem_categories) > 0) {
      filteredData <- data()[data()$Group %in% input$Chem_categories, ]
    } else {
      filteredData <- data()
    }
     
     filteredData
   })

  # Render the filtered results table
  output$resultsTable <- DT::renderDataTable({
    DT::datatable(filteredData(), 
                  options = list(
                    pageLength = 10,
                    columnDefs = list(
                      list(
                        targets = c(8,10:11, 17:21),
                        render = JS(
                          "function(data, type, row, meta) {",
                          "  if (type === 'display' || type === 'filter') {",
                          "    if (data === null || data === '') {",
                          "      return '';",
                          "    } else {",
                          "      if (meta.col >= 13 && meta.col <= 16) {",
                          "        return parseFloat(data).toExponential(2);",  # Convert to scientific notation with 2 decimals
                          "      } else {",
                          "        return parseFloat(data).toFixed(3);",  # Default behavior for other columns (3 decimals)
                          "      }",
                          "    }",
                          "  } else {",
                          "    return data;",
                          "  }",
                          "}")
                      ),
                      list(
                        targets = 2:3,
                        render = JS(
                          "function(data, type, row, meta) {",
                          "return type === 'display' && data.length > 20 ?",
                          "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                          "}"
                        )
                      )
                    )
                    ))
  }
  )
  
  # Clear the filter
  observeEvent(input$clearButton_1, {
    updateTextInput(session, "search", value = "")
  })
  
### Uncertainty-data ### -------------------------------------------------------
  # Filter the data based on user selections and display the results
  filteredData_uncert <- reactive({
    req(data_uncert())
    
    column_2 <- "CAS.Number"
    searchValue_2 <- input$search_2
    checkOnlyValid <- input$OnlyValid
    
    if (!is.null(column) && !is.null(searchValue_2)) {
      filteredData_uncert <- data_uncert()[grepl(searchValue_2, data_uncert()[[column_2]], ignore.case = TRUE), ]
    } else {
      filteredData_uncert <- data_uncert()
    }
    
    if(isTRUE(checkOnlyValid)){
     filteredData_uncert <- data_uncert() %>% filter(status == "Convergence")
    } else{
      filteredData_uncert <- data_uncert()
    }
    filteredData_uncert
  })
  
  # Render the filtered results table
  output$resultsTable_uncert <- renderDataTable({
    datatable(filteredData_uncert(), 
                  options = list(
                    pageLength = 25,
                    autoWidth = TRUE,
                    formatRound = 4,
                    dom = "ltirp",
                    columnDefs = list(
                      list(
                      targets = c(2:9, 12, 15),
                        render = JS(
                          "function(data, type, row, meta) {",
                          "  if (type === 'display' || type === 'filter') {",
                          "    if (data === null || data === '') {",
                          "      return '';",
                          "    } else {",
                          "      if ((meta.col >= 2 && meta.col <= 9) || meta.col === 15) {",
                          "        if (parseFloat(data) === 0) {",  # Check if value is zero
                          "          return '0';",  # Print as '0'
                          "        } else {",
                          "          return parseFloat(data).toExponential(2);",  # Convert to scientific notation with 2 decimals
                          "        }",
                          "      } else {",
                          "        if (meta.col === 12) {",
                          "          return parseFloat(data).toFixed(0);",  # Default behavior for other columns (3 decimals)
                          "        }",
                          "      }",
                          "    }",
                          "  } else {",
                          "    return data;",
                          "  }",
                          "}")
                      ))
                    ))
  }
  )
  
  # Clear the filter
  observeEvent(input$clearButton_2, {
    updateTextInput(session, "search_2", value = "")
  })
  
### SSD plots ### -------------------------------------------------
  # Read the NLS dataset
  valid_SSD_dataset <- read.csv("data/nls_output_df.csv") %>% filter(status == "Convergence") %>% pull(CAS.Number)
  
  SSDdata <- reactive({
 read.csv("data/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% 
      filter(CAS.Number %in% valid_SSD_dataset)
    })
  
  # # Update the CAS Number column choices
  # observeEvent(input$dataFile, {
  #   req(SSDdata(), input$casNumberInput, input$hcInput)
  #   updateSelectInput(session, "casNumberInput", choices = colnames(SSDdata()))
  # })
  
  # Run the function and generate the output
  outputData <- eventReactive(input$runButton, {
    req(SSDdata(), input$casNumberInput, input$hcInput, input$MC_n)
    nls_across_shiny(SSDdata(), CAS = input$casNumberInput, HCx = input$hcInput, MC_n = input$MC_n)
  })
  
  # Update the CAS Number choices based on the uploaded dataset
  observeEvent(SSDdata(), {
    updateSelectInput(session, "casNumberInput", choices = SSDdata() %>% distinct(CAS.Number) %>% pull(CAS.Number))
  })
  
  # Render the NLS output data table
  output$outputNLSData <- renderDataTable({
    # index the datatable object
    dat_list <- outputData()[[1]]
    outputDf <- data.frame(do.call(cbind, dat_list[1:16]), row.names = FALSE) %>%
        mutate(across(c(2:16), ~ as.numeric(.x))) %>% 
        `rownames<-`( NULL )
    return(outputDf)
  }, 
  options = list(
    dom = 'frtip'
    )
  )
  
  # Render the RAW data table used for SSDs
  output$outputRAWData <- renderDataTable({
    outputRAWDf <- SSDdata() %>% filter(CAS.Number == input$casNumberInput)
    return(outputRAWDf)
  },
  options = list(
    autoWidth = TRUE,
    #dom = 't',
    pageLength = 3000,
    columnDefs = list(
      list(
        targets = c(10,15:17),
        render = JS(
          "function(data, type, row, meta) {",
          "if (type === 'display' || type === 'filter') {",
          "  if (data === null || data === '') {",
          "    return '';",
          "  } else {",
          "    return parseFloat(data).toExponential(3);",  # Convert to scientific notation with 2 decimals
          "  }",
          "} else {",
          "  return data;",
          "}",
          "}")
      ),
      list(
        targets = 19,
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 10 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
          "}"
    ) ) ) )
  )
  
  # Render the output plot
  output$outputPlot <- renderPlot({
    req(input$casNumberInput, input$hcInput)
    # index the ggplot object
    outputData()[[2]]
    
  })
  
  # Render the output plot
  output$outputHist <- renderPlot({
    req(input$casNumberInput, input$hcInput, input$MC_n)
    # index the ggplot object
    hist_dat <- outputData()[[1]]
    hist(hist_dat$HCx_vec[[1]], main = paste("Histogram of Monte Carlo HC", hist_dat$Resp_lvl, "EC10eq-values. ", input$MC_n/1000, "K iterations", sep = ""))
    })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Tox_data_HESTIA_CASRN_", input$casNumberInput, ".csv", sep = "")
    },
    content = function(filename) {
      write.csv(SSDdata() %>% filter(CAS.Number == input$casNumberInput), 
                filename, 
                row.names = FALSE)
    }
  )
  
}

# Run the Shiny app
shinyApp(ui, server)
# , options = c(launch.browser = .rs.invokeShinyPaneViewer)