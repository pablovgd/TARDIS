library(shiny)
library(shinyFiles)
library(shinyjs)
library(jrc)
library(shinybusy)
library(bslib)
library(tidyverse)
library(readxl)
library(openxlsx)
library(MsExperiment)
library(BiocParallel)
library(Spectra)
library(RColorBrewer)
library(signal)
library(xcms)
library(pracma)


source(file = "R/peakdet.R")
source(file = "R/quality.R")
source(file = "R/createRanges.R")
source(file = "R/createTargetList.R")
source(file = "R/plotQCs.R")
source(file = "R/plotSamples.R")
source(file = "R/peaks_with_tardis.R")
source(file = "R/plotDiagnostic.R")
source(file = "R/sumintensities.R")

ui <- page_navbar(

  shinyjs::useShinyjs(),  # Initialize shinyjs
  title = div(img(src="/tardis.png",width =100)),
  bg = "#003B6F",
  inverse = TRUE,
  nav_panel(title = "Create target list",
            p(
              column(width = 6,
                     fileInput("target_file", "Table with target compounds (.csv or .xlsx):", accept = c(".xlsx",".csv")),
                     textInput("pos_pattern", "Positive pattern:", value = "+"),
                     textInput("neg_pattern", "Negative pattern:", value = "-"),
                     selectInput("polarity", "Polarity:", choices = c("positive","negative")),
                     textInput("ion_column", "Name of ion column:", value = "ion"),
                     textInput("columns_of_interest", "Names of columns of interest (comma-separated):", value = "id, name, mz, rt"),
                     actionButton("create_target_list", "Create Target List")
              )
            )),
  nav_panel(title = "Targeted peak detection",
            page_fillable(
              layout_columns(
                card(
                       shinyDirButton("dir", "Enter data folder path:",title ="directory "),
                       textOutput("selectedDir"),
                       shinyDirButton("dir_out", "Enter output folder path:",title ="directory out"),
                       textOutput("selectedDirOut")
                ),
                card(
                       numericInput("ppm", "PPM:", value = 5),
                       selectInput("mode", "Mode:", choices = c("metabolomics","lipidomics")),
                       uiOutput("mass_low"),
                       uiOutput("mass_high"),
                       textInput("rtdev", "RT deviation:", value = "18"),
                       textInput("int_std_id", "Internal standard IDs (comma-separated):", value = "331,1578,1576,1583,1577"),
                       textInput("batch_positions", "Batch postions (comma-separated):", value = "1,30"),
                       textInput("sample_pattern", "Sample pattern:", value = ""),
                       textInput("QC_pattern", "QC pattern:", value = "QC")
                       ),
                card(
                       checkboxInput("diagnostic_plots", "Generate diagnostic plots", value = TRUE),
                       checkboxInput("plot_samples", "Generate sample plots", value = FALSE),
                       checkboxInput("plot_QCs", "Generate QC plots", value = FALSE),
                       checkboxInput("screening_mode", "Screening mode", value = FALSE),
                       checkboxInput("batch_mode", "Batch mode", value = TRUE),
                       checkboxInput("rt_alignment", "RT alignment", value = TRUE),
                       checkboxInput("smoothing", "Smoothing", value = TRUE)
                       ),
                card(
                      actionButton("run_tardis_peaks", "Run T.A.R.D.I.S."),


                )
              )
            ))


    )

# Define server logic
server <- function(input, output, session) {

  shinyjs::useShinyjs()  # Initialize shinyjs

  volumes = getVolumes()()


  shinyDirChoose(input, "dir", roots = volumes)

  output$selectedDir <- renderText({
    if (!is.null(input$dir)) {
      paste("Selected Directory:", parseDirPath(volumes,input$dir))
    }
  })





  shinyDirChoose(input, "dir_out", roots = volumes)

  output$selectedDirOut <- renderText({
    if (!is.null(input$dir_out)) {
      paste("Selected Directory:", parseDirPath(volumes,input$dir_out))
    }
  })




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

    updateActionButton(inputId = "create_target_list", label = "Target list created")
    output$create_target_list_output <- renderPrint({
      targets$targets
    })
  })

  output$mass_low <- renderUI({
    if(input$mode == "lipidomics"){
      numericInput("mass_low","Mass range - lower limit",value = 67)
    } else {
      NULL
    }
  })

  output$mass_high <- renderUI({
    if(input$mode == "lipidomics"){
      numericInput("mass_high","Mass range - upper limit",value = 1000)
    } else {
      NULL
    }
  })



  # Function to run tardis_peaks
  observeEvent(input$run_tardis_peaks, {

    updateActionButton(inputId = "run_tardis_peaks", label = "Processing...")

    show_modal_spinner(
      spin = "double-bounce",
      color = "#003B6F",
      text = "Processing...",
      session = shiny::getDefaultReactiveDomain()
    )


    batch <- as.numeric(unlist(strsplit(input$batch_positions, ",")))

    result <- list()

    for (i in seq(1, length(batch), by = 2)) {
      result[[length(result) + 1]] <- c(batch[i], batch[i + 1])
    }


    tardis_output <- tardis_peaks(
      file_path = parseDirPath(volumes,input$dir),
      dbData = targets$targets,
      ppm = as.numeric(input$ppm),
      rtdev = as.numeric(input$rtdev),
      mode = input$mode,
      mass_range = c(input$mass_low,input$mass_high),
      screening_mode = input$screening_mode,
      polarity = input$polarity,
      output_directory = paste0(parseDirPath(volumes,input$dir_out),'/'),
      plots_samples = input$plot_samples,
      plots_QC = input$plot_QCs,
      smoothing = input$smoothing,
      diagnostic_plots = input$diagnostic_plots,
      batch_mode = input$batch_mode,
      batch_positions =  result,
      sample_pattern = input$sample_pattern,
      QC_pattern = input$QC_pattern,
      rt_alignment = input$rt_alignment,
      int_std_id = unlist(strsplit(input$int_std_id, ","))
    )


    remove_modal_spinner(session = getDefaultReactiveDomain())
    })
  }



# Run the application
shinyApp(ui = ui, server = server)
