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
           batch_mode,
           batch_positions,
           QC_pattern,
           sample_pattern) {
    
    results_samples <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0)
      )
    results_QCs <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0)
      )
    
    files <- list.files(file_path, full.names = T, pattern = "mzML")
    
    polarity <-
      polarity #you can change the polarity here, needs to be "positive" or "negative"
    mode <-
      mode #you can change the mode here, needs to be "lipidomics" or "metabolomics"
    #With these two parameters, you can set the allowed errors for mass & retention time!
    ppm <-
      ppm #ppm error for EIC extraction #POLAR = 5ppm / LIPIDOMICS =  10 ppm
    deltaTR = rtdev
    
    
    info_compounds <- dbData
    
    
    if (batch_mode == TRUE) {
      for (batchnr in 1:length(batch_positions)) {
        
        dbData <- info_compounds
        
        results_QCs_batch <-
          data.frame(
            Component = character(0),
            Sample = character(0),
            AUC = numeric(0),
            SNR = numeric(0),
            peak_cor = numeric(0),
            foundRT = numeric(0)
          )
        
        files_batch <-
          files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
        
        sample_files <-
          files_batch[grep(pattern = QC_pattern, files_batch, invert = TRUE)]
        QC_files <-
          files_batch[grep(pattern = QC_pattern, files_batch)]
        
        #First do the whole analysis for the QC's to find the RT time where the peaks are found.
        
        
        
        # start <- Sys.time()
        data_QC <-
          readMsExperiment(spectraFiles = QC_files, backend = MsBackendMzR(),BPPARAM =SnowParam(workers = 1))
        #data_QC <- Spectra(QC_files, source = MsBackendMzR())
        # stop <- Sys.time()
        # elapsed <- stop - start
        
        sample_names <-
          lapply(data_QC@sampleData$spectraOrigin, basename)
        
        ranges <- createRanges(data_QC, dbData, ppm, rtdev)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        
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
              filterRt(sample_spectra, rtRanges[j, ])
            sample_spectra <-
              filterMzRange(sample_spectra, mzRanges[j, ])
            
            sfs_agg <-
              addProcessing(sample_spectra, .sum_intensities)
            eic <-
              cbind(rtime(sfs_agg),
                    unlist(intensity(sfs_agg), use.names = FALSE))
            
            
            
            
            rt <- eic[, 1]
            int <- eic[, 2]
            
            int[which(is.na(int))] = 0
            
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            
            border <- find_peak_points(smoothed)
            
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
              
              found_rt <- rt[which(int == max(int))]
              max_int = max(int)
              
              # Get information about the current component from info_compounds
              compound_info <- dbData[j, ]
              
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
              compound_info <- dbData[j, ]
              
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
            plotQCs(compound_info,
                    output_directory,
                    rt_list,
                    int_list,
                    x_list,
                    y_list,
                    batchnr,
                    sample_names)
          }
          
        }
        
        results_QCs <- rbind(results_QCs, results_QCs_batch)
        #Next do the whole analysis for the samples in the same batch of the QC's to find the compounds at the corrected RT.
        
        
        data_samples <-
          readMsExperiment(spectraFiles = sample_files, backend = MsBackendMzR(),BPPARAM = SnowParam(workers=1))
        
        sample_names <-
          lapply(data_samples@sampleData$spectraOrigin, basename)
        
        
        # Replace rtmed with average foundRT from previous results
        
        new_rt_avg <- results_QCs_batch %>%
          group_by(ID) %>%
          summarise(mean = mean(foundRT), na.rm = TRUE)
        
        dbData <- merge(dbData,new_rt_avg,by = "ID")
        
        dbData$trold <- dbData$tr
        dbData$tr <- new_rt_avg$mean
        
        #If no RT is found, restore old RT
        
        for (k in 1:dim(dbData[1])) {
          if (is.na(dbData$tr[k]) == TRUE) {
            dbData$tr[k] <- dbData$trold[k]
          }
        }
        
        ranges <- createRanges(data_samples, dbData, ppm, rtdev)
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
              filterRt(sample_spectra, rtRanges[j, ])
            sample_spectra <-
              filterMzRange(sample_spectra, mzRanges[j, ])
            
            sfs_agg <-
              addProcessing(sample_spectra, .sum_intensities)
            eic <-
              cbind(rtime(sfs_agg),
                    unlist(intensity(sfs_agg), use.names = FALSE))
            
            
            
            
            rt <- eic[, 1]
            int <- eic[, 2]
            
            
            int[which(is.na(int))] = 0
            
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            
            border <- find_peak_points(smoothed)
            
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
              
              found_rt <- rt[which(int == max(int))]
              max_int = max(int)
              
              # Get information about the current component from info_compounds
              compound_info <- dbData[j, ]
              
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
              compound_info <- dbData[j, ]
              
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
            plotSamples(compound_info,
                        output_directory,
                        rt_list,
                        int_list,
                        x_list,
                        y_list,
                        batchnr,
                        sample_names)
          }
          
        }
      }
    }
    
    results <- results_samples
    
    auc_table <- results %>%
      select(Component, Sample, AUC) %>%
      spread(Sample, AUC)
    
    avg_metrics_table <- results %>%
      group_by(Component) %>%
      summarise_at(vars(-Sample), list(~ if (is.numeric(.))
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
    saveWorkbook(wb, paste0(output_directory, "feat_table.xlsx"))
    
    return(list(auc_table, avg_metrics_table))  
  }
