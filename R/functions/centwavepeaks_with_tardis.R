tardis_peaks_centwave <-
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
           sample_pattern,
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
    ppm <-
      ppm #ppm error for EIC extraction #POLAR = 5ppm / LIPIDOMICS =  10 ppm
    deltaTR = rtdev
    
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
     
      #Create ranges for all compounds
      ranges <- createRanges(data_QC, dbData, ppm, rtdev)
      #Get mz & rt ranges
      mzRanges <- ranges[[1]]
      rtRanges <- ranges[[2]]
      
      #Get the ranges for the internal standard compounds
      
      internal_standards_rt <- rtRanges[which(dbData$ID %in% int_std_id),]
      internal_standards_mz <- mzRanges[which(dbData$ID %in% int_std_id),]
      
      #Get QC sample names
      QC_sample_names <-
        lapply(data_QC@sampleData$spectraOrigin, basename)
      
      
      #Retrieve foundRT of internal standards in QC's using centWave with manualChromPeaks
      pks <- cbind(internal_standards_mz,internal_standards_rt)
      colnames(pks) <- c("mzmin","mzmax","rtmin","rtmax")
      res <- manualChromPeaks(data_QC,pks)
      chrompeaks_res <- data.frame(chromPeaks(res))
      incomplete <- which(table(chrompeaks_res$sample) < length(int_std_id) )
      chrompeaks_res <- chrompeaks_res[-which(chrompeaks_res$sample %in% incomplete),]
      rt_res <- chrompeaks_res$rt
      int_std <- matrix(rt_res,nrow = length(int_std_id),ncol=length(QC_sample_names)-length(incomplete),byrow = FALSE) 
      
      
      
      #Define parameters for retention time adjustment, based on the QC's and the internal standards
      param <- PeakGroupsParam(minFraction = 0.9, peakGroupsMatrix = int_std, subset =  which(sampleData(data_batch)$sample_type == "QC")[-incomplete])
      
      #Adjust the RT
      data_batch <- adjustRtime(data_batch, param = param)
      #Apply the adjusted RT 
      data_batch <- applyAdjustedRtime(data_batch)
      
      
      #Now, we try and find ALL compounds in the QC samples and save their foundRT to search the compounds at that RT in the sample files
      
      
      data_QC <- data_batch[which(sampleData(data_batch)$sample_type == "QC")]
      
      QC_sample_names <-
        lapply(data_QC@sampleData$spectraOrigin, basename)
      
      
   

      pks <- cbind(mzRanges,rtRanges)

      colnames(pks) <- c("mzmin", "mzmax", "rtmin" , "rtmax")

      res <- manualChromPeaks(data_QC, pks)
    
      #' pdp <- PeakDensityParam(sampleGroups = sampleData(res)$sample_type,
      #'                         minFraction = 0.1, bw = 10)
      #' res <- groupChromPeaks(res, param = pdp)
      #' 
      #' #' Get the feat peaks from our data set
      #' cpks <- as.data.frame(featureDefinitions(res))
      #' cpks$peak_id <- rownames(cpks)
      
      # matched <- matchValues(cpks,info_compounds,MzRtParam( ppm = 5, toleranceRt = 6),mzColname = c("mzmed","m/z"),rtColname = c("rtmed","tr"))
      # 
      # for(i in 1:dim(cpks)[1]){
      # cpks$comp_name[i] = matched@matches$target_idx[which(matched@matches$query_idx == i)]
      # }
      
      #for lusje voor plots
      dir.create(paste0(output_directory, "QCbatch_", batchnr))
      generate_QC_plot <- function(i) {
        res_chrom <- chromatogram(res, mz = mzRanges[i, ], rt = rtRanges[i, ],BPPARAM = bpparam())
        png(file = paste0(paste0(output_directory, "QCbatch_", batchnr, "/"),
                          paste(paste0("ID", info_compounds$ID[i]),".png")))
        plot(res_chrom, main = paste(paste0("ID", info_compounds$ID[i]), info_compounds$NAME[i]))
        dev.off()
        
      }
      
      # Use lapply to apply the function to each index
      future_lapply(1:length(info_compounds$NAME), generate_QC_plot)
      
      
      #Next do the whole analysis for the samples in the same batch of the QC's to find ALL the compounds at the corrected RT in both study samples & QC's.
      
      
   
      #Get sample data
      data_samples <- data_batch[which(sampleData(data_batch)$sample_type == "study")]
      
      
      sample_names <-
        lapply(data_samples@sampleData$spectraOrigin, basename)
      
      
      # Replace rtmed with average foundRT from previous results
      new_rt_avg <- c()
      chrompeaks <- as.data.frame(chromPeaks(res))
      
      for(i in 1:length(info_compounds$NAME)){ 
        pkscmp <- chrompeaks[which(chrompeaks$mzmin >= mzRanges[i,1] & chrompeaks$mzmax <= mzRanges[i,2] & chrompeaks$rtmin >= rtRanges[i,1] & chrompeaks$rtmax <= rtRanges[i,2]),]
        rt_avg <- mean(pkscmp$rt)
        new_rt_avg <- c(new_rt_avg, rt_avg)
      }
      
      
      
      dbData$trold <- dbData$tr
      dbData$tr <- new_rt_avg
      
      #If no RT is found, restore old RT
      
      for (k in 1:dim(dbData[1])) {
        if (is.na(dbData$tr[k]) == TRUE) {
          dbData$tr[k] <- dbData$trold[k]
        }
      }
      
      #Create ranges around new RT 
      ranges <- createRanges(data_samples, dbData, ppm, rtdev)
      mzRanges <- ranges[[1]]
      rtRanges <- ranges[[2]]
      
      
      
      pks <- cbind(mzRanges,rtRanges)
      
      colnames(pks) <- c("mzmin", "mzmax", "rtmin" , "rtmax")
      
      res <- manualChromPeaks(data_batch, pks)      
      
      feat_res <- groupChromPeaks(res,param = PeakDensityParam(sampleGroups = res@sampleData$sample_type))
      
      #Now we need to annotate the features to our target compounds
      
      matched <- matchValues(featureDefinitions(feat_res),dbData,MzRtParam( ppm = 5, toleranceRt = 6),mzColname = c("mzmed","m/z"),rtColname = c("rtmed","tr"))
      
      matched_feat_unique <- matched@matches[which(!duplicated(matched@matches$query_idx)),]
      
      feat_def_matched <- featureDefinitions(feat_res)[matched_feat_unique$query_idx,]
      
      feat_val_matched <- featureValues(feat_res)[matched_feat_unique$query_idx,]
      
      feat_def_matched$Name <- matched@target$NAME[matched_feat_unique$target_idx]
      feat_def_matched$ID <-  matched@target$ID[matched_feat_unique$target_idx]
      
      filtered_feat_res <- filterFeatureDefinitions(feat_res, features = rownames(feat_def_matched))
      
      chrom <- featureChromatograms(filtered_feat_res)
      
      for(i in 1:dim(chrom)[1]){
        plot(chrom[i,], main = paste(feat_def_matched$Name[i],"ID",feat_def_matched$ID[i]))
      }
      
      
    }
    
    results <- results_samples
    
    auc_table <- results %>%
      select(Component, Sample, AUC) %>%
      spread(Sample, AUC)
    
    avg_metrics_table <- results %>%
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
    saveWorkbook(wb, paste0(output_directory, "feat_table.xlsx"))
    
    return(list(auc_table, avg_metrics_table))
  }
