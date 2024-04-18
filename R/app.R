library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
source(file = "R/functions/peakdet.R")
source(file = "R/functions/quality.R")
source(file = "R/functions/createRanges.R")
source(file = "R/functions/createTargetList.R")
source(file = "R/functions/plotQCs.R")
source(file = "R/functions/plotSamples.R")
source(file = "R/functions/manualchrompeaks_fix.R")
source(file = "R/functions/peaks_with_tardis.R")
source(file = "R/functions/plotDiagnostic.R")
library(xcms)
library(Spectra)
library(MsExperiment)
library(RColorBrewer)
library(peakPantheR)
library(SummarizedExperiment)
library(MetaboAnnotation)
library(pracma)
library(signal)
library(tidyverse)
library(stringr)
library(readxl)
library(Polychrome)
library(openxlsx)
library(progress)
library(future)

ui <- fluidPage(
  titlePanel("T.A.R.D.I.S."),
  
  shinyjs::useShinyjs(),  # Initialize shinyjs
  
  sidebarLayout(
    sidebarPanel(
      # Inputs for the createTargetList function
      column(width = 6,
             fileInput("target_file", "Select Target File (xlsx format):", accept = c(".xlsx")),
             textInput("pos_pattern", "Positive Pattern:", value = "+"),
             textInput("neg_pattern", "Negative Pattern:", value = "-"),
             textInput("polarity", "Polarity:", value = "negative"),
             textInput("ion_column", "Ion Column:", value = "ion"),
             textInput("columns_of_interest", "Columns of Interest (comma-separated):", value = "id, name, mz, rt"),
             actionButton("create_target_list", "Create Target List")
      ),
      # Inputs for the tardis_peaks function
      column(width = 6,
             textInput("data_folder", "Enter Data Folder Path:"),
             textInput("ppm", "PPM:", value = "5"),
             textInput("rtdev", "RT Deviation:", value = "18"),
             textInput("output_directory", "Output Directory:", value = "D:/Data/PANIC_saliva/negative/"),
             checkboxInput("diagnostic_plots", "Generate Diagnostic Plots", value = TRUE),
             checkboxInput("batch_mode", "Batch Mode", value = TRUE),
             textInput("sample_pattern", "Sample Pattern:", value = ""),
             textInput("QC_pattern", "QC Pattern:", value = "QC"),
             checkboxInput("rt_alignment", "RT Alignment", value = TRUE),
             textInput("int_std_id", "Internal Standard IDs (comma-separated):", value = "1577, 1576, 1579"),
             actionButton("run_tardis_peaks", "Run Tardis Peaks"),
             div(id = "tardis_spinner", style = "display:none;", shinycssloaders::withSpinner(textOutput("tardis_processing_message")))
      )
    ),
    
    mainPanel(
      verbatimTextOutput("create_target_list_output"),
      verbatimTextOutput("tardis_peaks_output")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  shinyjs::useShinyjs()  # Initialize shinyjs
  
  # Reactive values to store targets
  targets <- reactiveValues()
  
  # Function to create target list
  observeEvent(input$create_target_list, {
    req(input$target_file)
    
    targets$targets <- createTargetList(
      input_directory_targets = input$target_file$datapath,
      pos_pattern = input$pos_pattern,
      neg_pattern = input$neg_pattern,
      polarity = input$polarity,
      ion_column = input$ion_column,
      columns_of_interest = unlist(strsplit(input$columns_of_interest, ", "))
    )
    
    output$create_target_list_output <- renderPrint({
      targets$targets
    })
  })
  
  # Function to run tardis_peaks
  observeEvent(input$run_tardis_peaks, {
    req(input$data_folder)
    
    # Show spinner while processing
    shinyjs::show("tardis_spinner")
    
    tardis_output <- tardis_peaks(
      file_path = input$data_folder,
      dbData = targets$targets,
      ppm = as.numeric(input$ppm),
      rtdev = as.numeric(input$rtdev),
      mode = "metabolomics",
      polarity = input$polarity,
      output_directory = input$output_directory,
      plots_samples = FALSE,
      plots_QC = FALSE,
      diagnostic_plots = input$diagnostic_plots,
      batch_mode = input$batch_mode,
      batch_positions = list(c(1, 208), c(209, 417)),
      sample_pattern = input$sample_pattern,
      QC_pattern = input$QC_pattern,
      rt_alignment = input$rt_alignment,
      int_std_id = as.numeric(unlist(strsplit(input$int_std_id, ", ")))
    )
    
    # Hide spinner after processing
    shinyjs::hide("tardis_spinner")
    
    output$tardis_peaks_output <- renderPrint({
      tardis_output
    })
  })
  
  # Reactive element for processing message
  output$tardis_processing_message <- renderText({
    "Processing..."
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
