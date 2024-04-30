tardis_peaks <-
  function(file_path,
           dbData,
           ppm,
           rtdev,
           mode,
           mass_range,
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
           int_std_id,
           screening_mode,
           smoothing) {
    #Create empty dataframe for sample results
    results_samples <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0),
        pop = numeric(0)
      )
    #Create empty dataframe for QC results
    results_QCs <-
      data.frame(
        Component = character(0),
        Sample = character(0),
        AUC = numeric(0),
        SNR = numeric(0),
        peak_cor = numeric(0),
        foundRT = numeric(0),
        pop = numeric(0)
      )
    
    #List files in directory
    files <-
      list.files(file_path, full.names = T, pattern = "mzML|mzXML")
    
    #Select polarity, mode, ppm & deltaTR
    polarity <-
      polarity #you can change the polarity here, needs to be "positive" or "negative"
    mode <-
      mode #you can change the mode here, needs to be "lipidomics" or "metabolomics"
    #With these two parameters, you can set the allowed errors for mass & retention time!
    ppm <-
      ppm #ppm error for EIC extraction #POLAR = 5ppm / LIPIDOMICS =  10 ppm
    
    mass_range <-
      mass_range
    
    deltaTR <- rtdev
    
    #save compound info
    if(mode == "lipidomics"){
      dbData <- dbData[which(dbData$`m/z` < mass_range[2] & dbData$`m/z` > mass_range[1]),]
    }
    
    
    info_compounds <- dbData
    
    
    
    
    if (screening_mode == TRUE) {
      #select all qc's
      QC_files <-
        files[grep(pattern = QC_pattern, files)]
      data_QC <- readMsExperiment(
        spectraFiles = QC_files,
        backend = MsBackendMzR(),
        BPPARAM = SnowParam(workers = 1)
      )
      
      if(mode == "lipidomics"){
        spectra_QC <- data_QC@spectra
        spectra_QC <- filterMzRange(spectra_QC,mass_range)
        spectra_QC <- filterEmptySpectra(spectra_QC)
      } else{spectra_QC <- data_QC@spectra}
      
      
      
      
      #Create ranges for all compounds
      ranges <- createRanges(data_QC, dbData, ppm, deltaTR)
      #Get mz & rt ranges
      mzRanges <- ranges[[1]]
      rtRanges <- ranges[[2]]
      
      
      if (rt_alignment == TRUE) {
        #Get the ranges for the internal standard compounds
        internal_standards_rt <-
          rtRanges[which(dbData$ID %in% int_std_id),]
        internal_standards_mz <-
          mzRanges[which(dbData$ID %in% int_std_id),]
        
        #Get QC sample names
        sample_names <-
          lapply(data_QC@sampleData$spectraOrigin, basename)
        
        
        #Initiate empty vectors
        int_std_foundrt <- c()
        int_std <- c()
        
        #Retrieve foundRT of internal standards in QC's, loop over all samples and all internal standards
        for (j in 1:dim(internal_standards_rt)[1]) {
          rt_list = list()
          int_list = list()
          x_list = list()
          y_list = list()
          
          
          
          for (i in 1:length(sample_names)) {
            #get current sample name
            sample_name <- unlist(sample_names[i])
            
            #Filter data to get data of sample and the mz and rt range of the current target
            sample_spectra <-
              filterDataOrigin(spectra_QC, unique(dataOrigin(spectra_QC))[i])
            sample_spectra <-
              filterRt(sample_spectra, internal_standards_rt[j,])
            sample_spectra <-
              filterMzRange(sample_spectra, internal_standards_mz[j,])
            
            #Extract intensities and RT from Spectra object
            sfs_agg <-
              addProcessing(sample_spectra, .sum_intensities)
            eic <-
              cbind(rtime(sfs_agg),
                    unlist(intensity(sfs_agg), use.names = FALSE))
            
            rt <- eic[, 1]
            int <- eic[, 2]
            
            #NA intensties are set to zero
            int[which(is.na(int))] = 0
            
            #To determine the borders more easily, smoothing is applied
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            
            
            if(smoothing == TRUE){
              int <- smoothed
              int[int < 0] <- 0
            }
            
            #Border detection
            border <- find_peak_points(rt, smoothed, dbData$tr[j])
            
            #Save found RT for internal standard target
            int_std_foundrt <-
              cbind(int_std_foundrt, border[[3]]) #this will finally contain all found rts for the different internal standards in this sample
            
          }
          
          
          int_std <-
            rbind(int_std, int_std_foundrt) #this will contain all the found rt's in all samples
          int_std_foundrt <- c()
        }
        
        
        #Define parameters for retention time adjustment, based on the QC's and the internal standards
        param <-
          PeakGroupsParam(
            minFraction = 0.9,
            span = 0.5  ,
            peakGroupsMatrix = int_std
          )
        
        #Adjust the RT
        data_QC <- adjustRtime(data_QC, param = param)
        #Apply the adjusted RT
        data_QC <- applyAdjustedRtime(data_QC)
      }
      
      #Find all targets in x QC's
      if(mode == "lipidomics"){
        spectra_QC <- data_QC@spectra
        spectra_QC <- filterMzRange(spectra_QC,mass_range)
        spectra_QC <- filterEmptySpectra(spectra_QC)
      } else{spectra_QC <- data_QC@spectra}
      
      spectra<- spectra_QC
      
      
      for (j in 1:dim(rtRanges)[1]) {
        rt_list = list()
        int_list = list()
        x_list = list()
        y_list = list()
        
        
        
        for (i in 1:length(sample_names)) {
          sample_name <- unlist(sample_names[i])
          
          
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
          
          
          if(smoothing == TRUE){
            int <- smoothed
            int[int < 0] <- 0
          }
          
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
            #calculate points over peak
            pop <- length(x)
            
            # Calculate QScore
            qscore <- qscoreCalculator(x, y)
            
            found_rt <- border$foundrt
            max_int = int[border$peakindex]
            
            # Get information about the current component from info_compounds
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_QCs <- rbind(
              results_QCs,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = auc,
                MaxInt = max_int,
                SNR = qscore[1],
                peak_cor = qscore[2],
                foundRT = found_rt,
                pop = pop,
                compound_info
              )
            )
            
          } else{
            compound_info <- dbData[j,]
            
            # Append results to the data frame with compound information
            results_QCs <- rbind(
              results_QCs,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = NA,
                MaxInt = NA,
                SNR = NA,
                peak_cor = NA,
                foundRT = NA,
                pop = NA,
                compound_info
              )
            )
          }
        }
        
        
        
        # Create and save the plot for the current component
        batchnr = 1
        if (diagnostic_plots == TRUE) {
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
      
      
      

      avg_metrics_table <- results_QCs %>%
        group_by(Component) %>%
        summarise_at(vars(-Sample), list( ~ if (is.numeric(.))
          mean(., na.rm = TRUE)
          else
            first(.)))
      
      write.csv(avg_metrics_table,
                file = paste0(output_directory, "qc_screening.csv"))
                
    } else {
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
            foundRT = numeric(0),
            pop = numeric(0)
          )
        
        #Select files in the batch
        files_batch <-
          files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
        
        #Subselection, only the sample files: use sample pattern if provided, if not, inverse of QC pattern
        
        if (sample_pattern == "") {
          sample_files <-
            files_batch[grep(pattern = QC_pattern, files_batch, invert = TRUE)]
        } else {
          sample_files <-
            files_batch[grep(pattern = sample_pattern, files_batch, invert = FALSE)]
        }
        #Subselection, only the QC files
        QC_files <-
          files_batch[grep(pattern = QC_pattern, files_batch)]
        
        #Load all the data for all files in the batch
        
        data_batch <- readMsExperiment(
          spectraFiles = files_batch,
          backend = MsBackendMzR(),
          BPPARAM = SnowParam(workers = 1)
        )
        
        #Define study and QC samples --> all not QC files are deemed study files
        sampleData(data_batch)$sample_type <- "study"
        sampleData(data_batch)$sample_type[grep(pattern = QC_pattern, files_batch)] <-
          "QC"
        
        #Subselect QC files, first we will locate the peaks of the internal standards in the QC files of this batch to align the rest of the batch data to
        
        data_QC <-
          data_batch[which(sampleData(data_batch)$sample_type == "QC")]
        #Extract spectra
        if(mode == "lipidomics"){
          spectra_QC <- data_QC@spectra
          spectra_QC <- filterMzRange(spectra_QC,mass_range)
          spectra_QC <- filterEmptySpectra(spectra_QC)
        } else{spectra_QC <- data_QC@spectra}
        #Create ranges for all compounds
        ranges <- createRanges(data_QC, dbData, ppm, deltaTR)
        #Get mz & rt ranges
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        
        
        if (rt_alignment == TRUE) {
          #Get the ranges for the internal standard compounds
          internal_standards_rt <-
            rtRanges[which(dbData$ID %in% int_std_id),]
          internal_standards_mz <-
            mzRanges[which(dbData$ID %in% int_std_id),]
          
          #Get QC sample names
          sample_names <-
            lapply(data_QC@sampleData$spectraOrigin, basename)
          
          
          #Initiate empty vectors
          int_std_foundrt <- c()
          int_std <- c()
          
          #Retrieve foundRT of internal standards in QC's, loop over all samples and all internal standards
          for (j in 1:dim(internal_standards_rt)[1]) {
            rt_list = list()
            int_list = list()
            x_list = list()
            y_list = list()
            
            
            
            for (i in 1:length(sample_names)) {
              #get current sample name
              sample_name <- unlist(sample_names[i])
              
              #Filter data to get data of sample and the mz and rt range of the current target
              sample_spectra <-
                filterDataOrigin(spectra_QC, unique(dataOrigin(spectra_QC))[i])
              sample_spectra <-
                filterRt(sample_spectra, internal_standards_rt[j,])
              sample_spectra <-
                filterMzRange(sample_spectra, internal_standards_mz[j,])
              
              #Extract intensities and RT from Spectra object
              sfs_agg <-
                addProcessing(sample_spectra, .sum_intensities)
              eic <-
                cbind(rtime(sfs_agg),
                      unlist(intensity(sfs_agg), use.names = FALSE))
              
              rt <- eic[, 1]
              int <- eic[, 2]
              
              #NA intensties are set to zero
              int[which(is.na(int))] = 0
              
              #To determine the borders more easily, smoothing is applied
              smoothed <- sgolayfilt(int, p = 3, n = 7)
              
              if(smoothing == TRUE){
                int <- smoothed
                int[int < 0] <- 0
              }
              
              #Border detection
              border <- find_peak_points(rt, smoothed, dbData$tr[j])
              
              #Save found RT for internal standard target
              int_std_foundrt <-
                cbind(int_std_foundrt, border[[3]]) #this will finally contain all found rts for the different internal standards in this sample
              
            }
            
            
            int_std <-
              rbind(int_std, int_std_foundrt) #this will contain all the found rt's in all samples
            int_std_foundrt <- c()
          }
          
          
          #Define parameters for retention time adjustment, based on the QC's and the internal standards
          param <-
            PeakGroupsParam(
              minFraction = 0.9,
              span = 0.5  ,
              peakGroupsMatrix = int_std,
              subset = which(sampleData(data_batch)$sample_type == "QC")
            )
          
          #Adjust the RT
          data_batch <- adjustRtime(data_batch, param = param)
          #Apply the adjusted RT
          data_batch <- applyAdjustedRtime(data_batch)
        }
        
        
        
        #Now, we try and find ALL compounds in the QC samples and save their foundRT to search the compounds at that RT in the sample files
        #Skip this step if there aren't any QC's available.
        
        
        data_QC <-
          data_batch[which(sampleData(data_batch)$sample_type == "QC")]
        
        if(length(data_QC) != 0 ){ 
        
          sample_names <-
            lapply(data_QC@sampleData$spectraOrigin, basename)
          
          
          
          if(mode == "lipidomics"){
            spectra_QC <- data_QC@spectra
            spectra_QC <- filterMzRange(spectra_QC,mass_range)
            spectra_QC <- filterEmptySpectra(spectra_QC)
          } else{spectra_QC <- data_QC@spectra}
          
          spectra <- spectra_QC
          
          for (j in 1:dim(rtRanges)[1]) {
            rt_list = list()
            int_list = list()
            x_list = list()
            y_list = list()
            
            
            
            for (i in 1:length(sample_names)) {
              sample_name <- unlist(sample_names[i])
              
              
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
              
              if(smoothing == TRUE){
                int <- smoothed
                int[int < 0] <- 0
              }
              
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
                
                pop <- length(x)
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
                    pop = pop,
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
                    pop = NA,
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
          
        }
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
        
        
      
        
        #Create ranges around new RT
        ranges <- createRanges(data_samples, dbData, ppm, deltaTR)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        
        
        if(mode == "lipidomics"){
          spectra <- data_samples@spectra
          spectra <- filterMzRange(spectra,mz = mass_range)
          spectra <- filterEmptySpectra(spectra)
        } else{spectra <- data_samples@spectra}

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
            
            
            eic <- eic[which(duplicated(eic[, 1]) == FALSE),]
            
            
            rt <- eic[, 1]
            int <- eic[, 2]
            
            
            int[which(is.na(int))] = 0
            
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            
            
            if(smoothing == TRUE){
              int <- smoothed
              int[int < 0] <- 0
            }
            
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
              
              pop <- length(x)
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
                  pop = pop,
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
                  pop = NA,
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
          
          if (diagnostic_plots == TRUE) {
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
      
      #AUC, int, SNR & peakcor tables for each component peak in every sample
      
      auc_table <- results %>%
        select(Component, Sample, AUC) %>%
        spread(Sample, AUC)
      
      write.csv(auc_table, file = paste0(output_directory, "auc_table.csv"))
      
      
      pop_table <- results %>%
        select(Component, Sample, pop) %>%
        spread(Sample, pop)
      
      write.csv(pop_table, file = paste0(output_directory, "pop_table.csv"))
      
      
      SNR_table <- results %>%
        select(Component, Sample, SNR) %>%
        spread(Sample, SNR)
      
      write.csv(SNR_table, file = paste0(output_directory, "snr_table.csv"))
      
      int_table <- results %>%
        select(Component, Sample, MaxInt) %>%
        spread(Sample, MaxInt)
      
      write.csv(int_table, file = paste0(output_directory, "int_table.csv"))
      
      
      peakcor_table <- results %>%
        select(Component, Sample, peak_cor) %>%
        spread(Sample, peak_cor)
      
      write.csv(peakcor_table, file = paste0(output_directory, "peakcor_table.csv"))
      
      
      #sumarrize feature table based on QC's
      
      avg_metrics_table <- NULL
      
      if(length(data_QC) != 0 ){ 
      
        QC_results <-  results[grep("QC", results$Sample),]
        
        avg_metrics_table <- QC_results %>%
          group_by(Component) %>%
          summarise_at(vars(-Sample), list( ~ if (is.numeric(.))
            mean(., na.rm = TRUE)
            else
              first(.)))
        
        
        
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
        saveWorkbook(wb,
                     paste0(output_directory, "feat_table.xlsx"),
                     overwrite = TRUE)
      }
      return(list(auc_table, avg_metrics_table))
    }
  }