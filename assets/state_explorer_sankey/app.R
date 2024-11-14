# libload
library(shiny)
library(plotly)
library(stringr)
library(readr)
library(shinyjs)

# Set max file upload size to 30 MB
options(shiny.maxRequestSize = 30 * 1024^2)

# Source custom Sankey plot function
source("sankey.R")

# User Interface
ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs for JavaScript functions
  
  # Hidden download button, to be triggered via JavaScript on click event
  conditionalPanel(
    "false", # Always hide this button in the UI
    downloadButton("downloadData")
  ),
  
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  titlePanel("Chromatin State Transitions"),
  
  # Sidebar layout with inputs and outputs
  sidebarLayout(
    sidebarPanel(
      
      # Input for uploading 'From' cell state file
      fileInput(
        inputId = "from_cell_bed",
        label = "From: Cell State",
        placeholder = "ChromHMM dense bed file",
        accept = ".bed"
      ),
      
      # Input for uploading 'To' cell state file
      fileInput(
        inputId = "to_cell_bed",
        label = "To: Cell State",
        placeholder = "ChromHMM dense bed file",
        accept = ".bed"
      ),
      
      tags$hr(),
      
      # Dynamic dropdown for selecting the state transition
      uiOutput("fromstate")
    ),
    
    mainPanel(
      
      # Output for information and error messages
      htmlOutput("info_handler"),
      htmlOutput("error_handler"),
      
      # Output for the Sankey plot
      plotlyOutput("sankey_plot"),
      
      # Output for displaying details of clicked plot elements
      verbatimTextOutput("click")
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive function to read 'From' cell state file data
  filedata <- reactive({
    infile <- input$from_cell_bed
    if (is.null(infile)) return(NULL)
    
    # Read the input file and set column names
    read_tsv(infile$datapath, skip = 1, col_names = c("chr", "from", "to", "statenumber", "score", "start", "thickon", "thickoff", "itemrgb"))
  })
  
  # Dynamic UI to populate 'State for transition' dropdown based on input data
  output$fromstate <- renderUI({
    df <- filedata()
    if (is.null(df)) return(NULL)
    
    items <- sort(unique(df$statenumber))
    names(items) <- items
    selectInput("fromstate", "State for transition:", items)
  })
  
  # Render the Sankey plot based on selected inputs
  output$sankey_plot <- renderPlotly({
    validate(
      need(input$from_cell_bed, ""),
      need(input$to_cell_bed, ""),
      need(input$fromstate, "")
    )
    
    # Generate Sankey plot using custom plot_sankey function
    plot_sankey(
      from_cell_bed = input$from_cell_bed,
      to_cell_bed = input$to_cell_bed,
      fromstate = input$fromstate
    )
  })
  
  # Event handler for clicks on the Sankey plot
  observeEvent(event_data("plotly_click"), {
    barData <- event_data("plotly_click")
    runjs("$('#downloadData')[0].click();")  # Trigger download on click
    to_state_clicked_id <<- barData[['pointNumber']] + 1  # Capture clicked state ID
  })
  
  # Download handler to save specific transition data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("regions-from-", input$fromstate, "-to-", to_state_clicked_id, ".bed", sep = "")
    },
    content = function(file) {
      # System command to filter and save data based on selected states
      command_run_outbed <- paste0(
        'tail -n +2 ', input$from_cell_bed[['datapath']],
        ' | cut -f 1,2,3,4 | ggrep -P "\\t', input$fromstate,
        '$" | bedtools intersect -a stdin -b ', input$to_cell_bed[['datapath']],
        ' -wao | ggrep -v -P "\\t\\.\\t\\d+$" | ggrep -P "\\t', to_state_clicked_id,
        '\\t" | bedtools sort | uniq > ', file
      )
      system(command_run_outbed)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)