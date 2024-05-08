server <- function(input, output, session) {

  shinyjs::useShinyjs()  # Initialize shinyjs

  volumes = shinyFiles::getVolumes()()


  shinyFiles::shinyDirChoose(input, "dir", roots = volumes)

  output$selectedDir <- renderText({
    if (!is.null(input$dir)) {
      paste("Selected Directory:", shinyFiles::parseDirPath(volumes,input$dir))
    }
  })





  shinyFiles::shinyDirChoose(input, "dir_out", roots = volumes)

  output$selectedDirOut <- renderText({
    if (!is.null(input$dir_out)) {
      paste("Selected Directory:", shinyFiles::parseDirPath(volumes,input$dir_out))
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

    updateActionButton(inputId = "run_tardis_peaks", label = "Processing done! You may close this window.")

    shinybusy::show_modal_spinner(
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
      file_path = shinyFiles::parseDirPath(volumes,input$dir),
      dbData = targets$targets,
      ppm = as.numeric(input$ppm),
      rtdev = as.numeric(input$rtdev),
      mode = input$mode,
      mass_range = c(input$mass_low,input$mass_high),
      screening_mode = input$screening_mode,
      polarity = input$polarity,
      output_directory = paste0(shinyFiles::parseDirPath(volumes,input$dir_out),'/'),
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


    shinybusy::remove_modal_spinner(session = getDefaultReactiveDomain())
  })
}
