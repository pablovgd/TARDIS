tardis_peaks <-
  function(file_path,
           dbData,
           ppm,
           rtdev,
           mode,
           polarity,
           output_directory,
           plots_samples,
           plots_QC,
           diagnostic_plots,
           batch_mode,
           batch_positions,
           QC_pattern,
           sample_pattern,
           rt_alignment,
           int_std_id) {
    #Create empty dataframe for sample results
    results_samples <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0)
      )
    #Create empty dataframe for QC results
    results_QCs <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0)
      )
    
    #List files in directory
    files <- list.files(file_path, full.names = T, pattern = "mzML")
    
    #Select polarity, mode, ppm & deltaTR
    polarity <-
      polarity #you can change the polarity here, needs to be "positive" or "negative"
    mode <-
      mode #you can change the mode here, needs to be "lipidomics" or "metabolomics"
    #With these two parameters, you can set the allowed errors for mass & retention time!
    ppm <- ppm #ppm error for EIC extraction #POLAR = 5ppm / LIPIDOMICS =  10 ppm
    
    deltaTR <- rtdev
    
    #save compound info
    info_compounds <- dbData
    
    #Loop over the batches
    
    for (batchnr in 1:length(batch_positions)) {
      dbData <- info_compounds
      
      #Empty dataframe for results of the QCs
      results_QCs_batch <-
        data.frame(
          Component = character(0),
          Sample = character(0),
          AUC = numeric(0),
          SNR = numeric(0),
          peak_cor = numeric(0),
          foundRT = numeric(0)
        )
      
      #Select files in the batch
      files_batch <-
        files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
      
      #Subselection, only the sample files
      sample_files <-
        files_batch[grep(pattern = QC_pattern, files_batch, invert = TRUE)]
      #Subselection, only the sample files 
      QC_files <-
        files_batch[grep(pattern = QC_pattern, files_batch)]
      
      #Load all the data for all files in the batch
      
      data_batch <- readMsExperiment(spectraFiles = files_batch,
                                     backend = MsBackendMzR(),
                                     BPPARAM = SnowParam(workers = 1)
      )
      
      #Define study and QC samples
      sampleData(data_batch)$sample_type <- "study"
      sampleData(data_batch)$sample_type[grep(pattern = QC_pattern, files_batch)] <- "QC"
      
      #Subselect QC files, first we will locate the peaks of the internal standards in the QC files of this batch to align the rest of the batch data to
      
      data_QC <- data_batch[which(sampleData(data_batch)$sample_type == "QC")]
      #Extract spectra
      spectra_QC <- data_QC@spectra
      #Create ranges for all compounds
      ranges <- createRanges(data_QC, dbData, ppm, deltaTR)
      #Get mz & rt ranges
      mzRanges <- ranges[[1]]
      rtRanges <- ranges[[2]]
      
      
      if(rt_alignment == TRUE){
        #Get the ranges for the internal standard compounds
        
        
        internal_standards_rt <- rtRanges[which(dbData$ID %in% int_std_id),]
        internal_standards_mz <- mzRanges[which(dbData$ID %in% int_std_id),]
        
        #Get QC sample names
        sample_names <-
          lapply(data_QC@sampleData$spectraOrigin, basename)
        
        
        #Initiate empty vectors
        int_std_foundrt <- c()
        int_std <- c()
        
        #Retrieve foundRT of internal standards in QC's
        for (j in 1:dim(internal_standards_rt)[1]) {
          rt_list = list()
          int_list = list()
          x_list = list()
          y_list = list()
          
          
          
          for (i in 1:length(sample_names)) {
            sample_name <- unlist(sample_names[i])
            
            #extraction of rt & int with chromatograms (slow)
            # chromatograms <-
            #   chromatogram(data_QC, rt = rtRanges[j, ], mz = mzRanges[j, ],aggregationFun = "max")
            # rt <- chromatograms@.Data[[i]]@rtime
            # int <- chromatograms@.Data[[i]]@intensity
            
            #extraction of rt & int using spectra (fast, but diff results??)
            
            
            .sum_intensities <- function(x, ...) {
              if (nrow(x)) {
                cbind(mz = NA_real_,
                      intensity = sum(x[, "intensity"], na.rm = TRUE))
              } else
                cbind(mz = NA_real_, intensity = NA_real_)
            }
            
            sample_spectra <-
              filterDataOrigin(spectra_QC, unique(dataOrigin(spectra_QC))[i])
            sample_spectra <-
              filterRt(sample_spectra, internal_standards_rt[j,])
            sample_spectra <-
              filterMzRange(sample_spectra, internal_standards_mz[j,])
            
            sfs_agg <-
              addProcessing(sample_spectra, .sum_intensities)
            eic <-
              cbind(rtime(sfs_agg),
                    unlist(intensity(sfs_agg), use.names = FALSE))
            
            
            
            
            rt <- eic[, 1]
            int <- eic[, 2]
            
            int[which(is.na(int))] = 0
            
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            
            border <- find_peak_points(rt, smoothed,dbData$tr[j])
            
            int_std_foundrt <- cbind(int_std_foundrt, border[[3]])
            
          }
          
          int_std <- rbind(int_std,int_std_foundrt)
          int_std_foundrt <- c()
        }
        
        
        # rtmed <- dbData$tr
        # param <- ObiwarpParam(binSize = 0.1,subset = which(sampleData(data_batch)$sample_type == "QC"))
        # rts <- as.matrix(cbind(rtmed,rtmed,rtmed,rtmed))
        # internal_standards <- rts[c(19,32,30,33,31),]
        
        #Define parameters for retention time adjustment, based on the QC's and the internal standards
        param <- PeakGroupsParam(minFraction = 0.9,span = 0.5  ,peakGroupsMatrix = int_std, subset = which(sampleData(data_batch)$sample_type == "QC"))
        
        #Adjust the RT
        data_batch <- adjustRtime(data_batch, param = param)
        #Apply the adjusted RT 
        data_batch <- applyAdjustedRtime(data_batch)
      }
      
      #Now, we try and find ALL compounds in the QC samples and save their foundRT to search the compounds at that RT in the sample files
      
      
      data_QC <- data_batch[which(sampleData(data_batch)$sample_type == "QC")]
      
      sample_names <-
        lapply(data_QC@sampleData$spectraOrigin, basename)
      
    
      
      spectra <- data_QC@spectra
 
      
      for (j in 1:dim(rtRanges)[1]) {
        rt_list = list()
        int_list = list()
        x_list = list()
        y_list = list()
        
        
        
        for (i in 1:length(sample_names)) {
          sample_name <- unlist(sample_names[i])
          
          #extraction of rt & int with chromatograms (slow)
          # chromatograms <-
          #   chromatogram(data_QC, rt = rtRanges[j, ], mz = mzRanges[j, ],aggregationFun = "max")
          # rt <- chromatograms@.Data[[i]]@rtime
          # int <- chromatograms@.Data[[i]]@intensity
          
          #extraction of rt & int using spectra (fast, but diff results??)
          
          
          .sum_intensities <- function(x, ...) {
            if (nrow(x)) {
              cbind(mz = NA_real_,
                    intensity = sum(x[, "intensity"], na.rm = TRUE))
            } else
              cbind(mz = NA_real_, intensity = NA_real_)
          }
          
          sample_spectra <-
            filterDataOrigin(spectra, unique(dataOrigin(spectra))[i])
          sample_spectra <-
            filterRt(sample_spectra, rtRanges[j,])
          sample_spectra <-
            filterMzRange(sample_spectra, mzRanges[j,])
          
          sfs_agg <-
            addProcessing(sample_spectra, .sum_intensities)
          eic <-
            cbind(rtime(sfs_agg),
                  unlist(intensity(sfs_agg), use.names = FALSE))
          
          
          
          
          rt <- eic[, 1]
          int <- eic[, 2]
          
          int[which(is.na(int))] = 0
          
          smoothed <- sgolayfilt(int, p = 3, n = 7)
          
          border <- find_peak_points(rt, smoothed,dbData$tr[j])
          
          x <- rt[border$left:border$right]
          y <- int[border$left:border$right]
          
          rt_list <- c(rt_list, list(rt))
          int_list <- c(int_list, list(int))
          x_list <- c(x_list, list(x))
          y_list <- c(y_list, list(y))
          
          # Check if there are at least two unique values for the component
          if (length(unique(y)) > 1) {
            # Calculate AUC
            auc <- trapz(x, y)
            
            # Calculate QScore
            qscore <- qscoreCalculator(x, y)
            
            found_rt <- border$foundrt
            max_int = int[border$peakindex]
            
            # Get information about the current component from info_compounds
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_QCs_batch <- rbind(
              results_QCs_batch,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = auc,
                MaxInt = max_int,
                SNR = qscore[1],
                peak_cor = qscore[2],
                foundRT = found_rt,
                compound_info
              )
            )
            
          } else{
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_QCs_batch <- rbind(
              results_QCs_batch,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = NA,
                MaxInt = NA,
                SNR = NA,
                peak_cor = NA,
                foundRT = NA,
                compound_info
              )
            )
          }
        }
        
        
        
        # Create and save the plot for the current component
        
        if (plots_QC == TRUE) {
          plotQCs(
            compound_info,
            output_directory,
            rt_list,
            int_list,
            x_list,
            y_list,
            batchnr,
            sample_names
          )
        }
        
      }
      
      results_QCs <- rbind(results_QCs, results_QCs_batch)
      
      #Next do the whole analysis for the samples in the same batch of the QC's to find ALL the compounds at the corrected RT. (SAMPLES + QC)
      
      
      # data_samples <-
      #   readMsExperiment(
      #     spectraFiles = sample_files,
      #     backend = MsBackendMzR(),
      #     BPPARAM = SnowParam(workers = 1)
      #   )
      
      #Get sample data
      data_samples <- data_batch
      #data_samples <- data_batch[which(sampleData(data_batch)$sample_type == "study")]
      
      
      # data_samples <- adjustRtime(data_samples,param = ObiwarpParam()) 
    
      
      sample_names <-
        lapply(data_samples@sampleData$spectraOrigin, basename)
      
      
      # Replace rtmed with average foundRT from previous results
      
      new_rt_avg <- results_QCs_batch %>%
        group_by(ID) %>%
        summarise(mean = mean(foundRT), na.rm = TRUE)
      
      dbData <- merge(dbData, new_rt_avg, by = "ID")
      
      dbData$trold <- dbData$tr
      dbData$tr <- new_rt_avg$mean
      
      #If no RT is found, restore old RT
      
      for (k in 1:dim(dbData[1])) {
        if (is.na(dbData$tr[k]) == TRUE) {
          dbData$tr[k] <- dbData$trold[k]
        }
      }
      
      #Create ranges around new RT 
      ranges <- createRanges(data_samples, dbData, ppm, deltaTR)
      mzRanges <- ranges[[1]]
      rtRanges <- ranges[[2]]
      
      spectra <- data_samples@spectra
      
      for (j in 1:dim(rtRanges)[1]) {
        rt_list = list()
        int_list = list()
        x_list = list()
        y_list = list()
        
        for (i in 1:length(sample_names)) {
          sample_name <- unlist(sample_names[i])
          # chromatograms <-
          #   chromatogram(data_QC, rt = rtRanges[j, ], mz = mzRanges[j, ])
          # rt <- chromatograms@.Data[[i]]@rtime
          # int <- chromatograms@.Data[[i]]@intensity
          
          
          sample_spectra <-
            filterDataOrigin(spectra, unique(dataOrigin(spectra))[i])
          sample_spectra <-
            filterRt(sample_spectra, rtRanges[j,])
          sample_spectra <-
            filterMzRange(sample_spectra, mzRanges[j,])
          
          sfs_agg <-
            addProcessing(sample_spectra, .sum_intensities)
          eic <-
            cbind(rtime(sfs_agg),
                  unlist(intensity(sfs_agg), use.names = FALSE))
          
          
          eic <- eic[which(duplicated(eic[,1]) == FALSE),]
          
          
          rt <- eic[, 1]
          int <- eic[, 2]
          
          
          int[which(is.na(int))] = 0
          
          smoothed <- sgolayfilt(int, p = 3, n = 7)
          
          border <- find_peak_points(rt, smoothed, dbData$tr[j])
          
          x <- rt[border$left:border$right]
          y <- int[border$left:border$right]
          
          rt_list <- c(rt_list, list(rt))
          int_list <- c(int_list, list(int))
          x_list <- c(x_list, list(x))
          y_list <- c(y_list, list(y))
          
          # Check if there are at least two unique values for the component
          if (length(unique(y)) > 1) {
            # Calculate AUC
            auc <- trapz(x, y)
            
            # Calculate QScore
            qscore <- qscoreCalculator(x, y)
            
            found_rt <- border$foundrt
            max_int = int[border$peakindex]
            
            # Get information about the current component from info_compounds
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_samples <- rbind(
              results_samples,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = auc,
                MaxInt = max_int,
                SNR = qscore[1],
                peak_cor = qscore[2],
                foundRT = found_rt,
                compound_info
              )
            )
            
          } else{
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_samples <- rbind(
              results_samples,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = NA,
                MaxInt = NA,
                SNR = NA,
                peak_cor = NA,
                foundRT = NA,
                compound_info
              )
            )
          }
          
        }
        # Create and save the plot for the current component
        
        if (plots_samples == TRUE) {
          plotSamples(
            compound_info,
            output_directory,
            rt_list,
            int_list,
            x_list,
            y_list,
            batchnr,
            sample_names
          )
        }
        
        if (diagnostic_plots == TRUE){
          plotDiagnostic(
            compound_info,
            output_directory,
            rt_list,
            int_list,
            x_list,
            y_list,
            batchnr,
            sample_names
          )
        }
      }
    }
    
    
    results <- results_samples
    
    auc_table <- results %>%
      select(Component, Sample, AUC) %>%
      spread(Sample, AUC)
    
    #sumarrize feature table based on QC's
    
    QC_results <-  results[grep("QC",results$Sample),]
    
    avg_metrics_table <- QC_results %>%
      group_by(Component) %>%
      summarise_at(vars(-Sample), list( ~ if (is.numeric(.))
        mean(., na.rm = TRUE)
        else
          first(.)))
    
    
    write.csv(auc_table, file = paste0(output_directory, "auc_table.csv"))
    
    lower_threshold <- 0.7
    upper_threshold <- 0.8
    
    # Create a new Excel workbook
    wb <- createWorkbook()
    
    # Add a worksheet to the workbook
    sheet <- addWorksheet(wb, "Sheet1")
    
    # Write the data.frame to the Excel worksheet
    writeData(wb,
              sheet,
              avg_metrics_table,
              startCol = 1,
              startRow = 1)
    
    # Define colors for formatting
    color_low <- "#FF7F7F"
    color_high <- "#90EE90"
    color_mid <- "#FFD580"
    
    # Apply conditional formatting based on the values in the "Score" column
    for (i in 1:nrow(avg_metrics_table)) {
      score <- avg_metrics_table$peak_cor[i]
      if (is.na(score) == TRUE) {
        score = 0
      }
      
      if (score < lower_threshold) {
        addStyle(wb,
                 sheet,
                 createStyle(fgFill = color_low) ,
                 i + 1 ,
                 1:ncol(avg_metrics_table))
      } else if (score > upper_threshold) {
        addStyle(wb,
                 sheet,
                 createStyle(fgFill = color_high) ,
                 i + 1 ,
                 1:ncol(avg_metrics_table))
      } else{
        addStyle(wb,
                 sheet,
                 createStyle(fgFill = color_mid) ,
                 i + 1 ,
                 1:ncol(avg_metrics_table))
      }
    }
    
    # Save the Excel workbook to a file
    saveWorkbook(wb, paste0(output_directory, "feat_table.xlsx"),overwrite = TRUE)
    
    return(list(auc_table, avg_metrics_table))
  }
