#' Find targeted compounds in LC-MS data
#'
#' @param path_data_in path to the .mzML or .mzXML files
#' @param msRawData_samples output of selectSpectra
#' @param msRawData_QC output of selectSpectra
#' @param QC_list list indicating the QC files for each batch
#' @param batch_list list indicating the different batches
#' @param ppm Allowed mass deviation in ppm
#' @param deltaTR Allowed RT deviation in the search window
#' @param dbData Output of loadTargetedDatabase
#'
#' @return Returns intensities/areas for found targeted compounds.
#' @export
#'
#' @examples
findTargetedCompounds <-
  function(path_data_in,
           msRawData_samples,
           msRawData_QC,
           QC_list,
           batch_list,
           ppm,
           deltaTR,
           dbData,
           plots) {
    trRange = 30
    m = 200
    
    
    
    
    
    #RUN
    for (counter in 1:length(QC_list)) {
      setwd(path_data_in)
      files <- QC_list[[counter]]
      
      #LOAD RAW DATA OF QC's OF CHOICE --> Choose to QC's you want to include
      
      xcmsset <- xcmsSet(msRawData_QC[files])
      msRawData_QCs <-
        getXcmsRaw(xcmsset, sampleidx = seq(1, length(files), 1))
      
      
      
      mergedata_list_QC <- list()
      targExtracRes_list_QC <- list()
      
      #RUN THIS LOOP
      for (p in 1:length(msRawData_QCs)) {
        #Load raw data
        rawData <- msRawData_QCs[[p]]
        
        #select mz & rt from database
        mzmed <- dbData$`m/z`
        rtmed <- dbData$tr
        #combine m/z & rt & give unique identifier
        mzAndRt <-
          as.data.frame(cbind(c(1:length(mzmed)), paste(mzmed, rtmed)))
        colnames(mzAndRt) <- c("ID", "mzAndRt")
        mzAndRt$ID <- as.numeric(mzAndRt$ID)
        
        #find unique combinations of mz & RT and give identifier
        uniqueMzAndRtID <- which(!duplicated(mzAndRt$mzAndRt))
        uniqueMzAndRt <- mzAndRt[uniqueMzAndRtID, ]
        uniqueMzAndRt$ID <- c(1:nrow(uniqueMzAndRt))
        colnames(uniqueMzAndRt) <- c("uniqueID", "mzAndRt")
        
        #merge mz & RT combinations with their unique ID
        mergeData <- merge(mzAndRt, uniqueMzAndRt, by = "mzAndRt")
        mergeData_QC <- mergeData[order(mergeData$ID), ]
        
        
        #save mz & rt from the unique mz & RT combo's
        mzmed <- mzmed[uniqueMzAndRtID]
        rtmed <- rtmed[uniqueMzAndRtID]
        
        
        #calculate delta m/z for each mass based on input ppm
        mzdeltas <- sapply(mzmed, function(mzmed)
          mzmed * ppm / 10 ^ 6)
        
        #calculate mzrange based on delta mz
        mzRanges <-
          cbind(as.numeric(mzmed) - mzdeltas,
                as.numeric(mzmed) + mzdeltas)
        
        #if an upper m/z boundary is lower than the minimum m/z range, set it to the minimum m/z
        indexTemp <- which(mzRanges[, 2] < min(rawData@mzrange))
        mzRanges[indexTemp, ] <- min(rawData@mzrange)
        
        #if a lower m/z boundary is higher than the max m/z, set it to the max m/z
        indexTemp <- which(mzRanges[, 1] > max(rawData@mzrange))
        mzRanges[indexTemp, ] <- max(rawData@mzrange)
        
        #if an upper limit is higher than the max mz & the lower limit is smaller than the max m/z, set the upper limit to the max m/z
        mzRanges[which(mzRanges[, 2] > max(rawData@mzrange) &
                         mzRanges[, 1] < max(rawData@mzrange)), 2] <- max(rawData@mzrange)
        
        #if a upper limit is larger than the minimum & the lower limit is lower than the minimum, set the lower limit to the min m/z
        mzRanges[which(mzRanges[, 2] > min(rawData@mzrange) &
                         mzRanges[, 1] < min(rawData@mzrange)), 1] <- min(rawData@mzrange)
        
        
        #calculate rt range based on deltaTR
        rtRanges <-
          cbind(as.numeric(rtmed) - deltaTR / 2,
                as.numeric(rtmed) + deltaTR / 2)
        
        #for al upper RT limits -2 lower than minimum RT, change to lower limit of min rt & upper limit to min + 10
        indexTemp <- which(rtRanges[, 2] - 2 < min(rawData@scantime))
        rtRanges[indexTemp, ] <-
          cbind(rep(min(rawData@scantime), length(indexTemp)), rep(min(rawData@scantime) +
                                                                     10, length(indexTemp)))
        
        
        #for al lower RT limits +2 higher than maximum RT, change lower limit to max -10 and upper limit to max
        indexTemp <- which(rtRanges[, 1] + 2 > max(rawData@scantime))
        rtRanges[indexTemp, ] <-
          cbind(rep(max(rawData@scantime) - 10, length(indexTemp)), rep(max(rawData@scantime), length(indexTemp)))
        
        
        #for upper limits higher than max and lower limit lower than max, change upper limit to max
        rtRanges[which(rtRanges[, 2] > max(rawData@scantime) &
                         rtRanges[, 1] < max(rawData@scantime)), 2] <- max(rawData@scantime)
        #for upper limits higher than min and lower limit lower than min, change lower limit to min
        rtRanges[which(rtRanges[, 2] > min(rawData@scantime) &
                         rtRanges[, 1] < min(rawData@scantime)), 1] <- min(rawData@scantime)
        
        #extract EIC data for defined ranges
        # EICdata <- getEIC(rawData, mzrange = mzRanges, rtrange = rtRanges) #this function uses profile matrix data
        
        rawEICdata_list <- list()
        for (i in 1:dim(rtRanges)[1]) {
          rawEICdata <-
            plotEIC(rawData, mzrange = mzRanges[i, ], rtrange = rtRanges[i, ]) #using this we can extract raw data
          rawEICdata_list[[i]] <- rawEICdata
        }
        
        
        
        #extract raw data
        # EICdataEIC <- EICdata@eic$xcmsRaw
        # #select the ones that have a uniqueID
        # EICdataEIC <- EICdataEIC[mergeData$uniqueID]
        rawEICdata_list <- rawEICdata_list[mergeData_QC$uniqueID]
        
        
        
        
        
        #list to char function to create vectors with all EIC data per mzrange & rtrange
        list2Character <- function(ithMatrix) {
          return(paste(paste(ithMatrix[, 1], ithMatrix[, 2], sep = " "), collapse = ";"))
        }
        
        
        
        # EICdataEICVector <- unlist(lapply(EICdataEIC, list2Character))
        rawEICdatavector <-
          unlist(lapply(rawEICdata_list, list2Character))
        #add to dbdata
        dbData$EIC <- rawEICdatavector
        
        #functie maken die "peak detection" en entrophy calc toepast op elke set EIC.
        # binnen elke EIC lokale maxima selecteren en als output RT & hoogte van lokaal maximum
        #er zit geen functie voor area te berekenen in de functie
        #berekend entropie binnen trrange interval
        func <- function(ithRowDbData, ithEICdataEIC) {
          extractedPeaks <-
            peakDectAndEntroCal(ithEICdataEIC, trRange = trRange, m = m)
          return(cbind(ithRowDbData[rep(1, nrow(extractedPeaks)), ], extractedPeaks))
        }
        
        dbDataList <- split(dbData, 1:nrow(dbData))
        
        targExtracRes_QC <-
          mapply(func, dbDataList, rawEICdata_list, SIMPLIFY = F)
        targExtracRes_QC <-
          do.call(rbind, lapply(targExtracRes_QC, data.frame))
        targExtracRes_QC$trOfPeak <-
          as.numeric(targExtracRes_QC$trOfPeak)
        
        
        
        mergedata_list_QC <- c(mergedata_list_QC, list(mergeData_QC))
        targExtracRes_list_QC <-
          c(targExtracRes_list_QC, list(targExtracRes_QC))
        
        
        
        
        
        
      }
      
      
      
      
      int_table_QC <- data.frame(targExtracRes_list_QC[[1]]$NAME)
      
      
      for (k in 1:length(targExtracRes_list_QC)) {
        int_table_QC <-
          cbind(int_table_QC, targExtracRes_list_QC[[k]]$peakHeight)
      }
      
      colnames(int_table_QC) <- append("Name", msRawData_QC[files])
      
      
      
      #PLOTS for the QC
      
      #OPTIONAL IF YOU WANT TO HAVE THE PLOTS IN A DIFFERENT DIR YOU CAN CHANGE IT HERE
      dir.create(paste0(path_data_in, "/batch_", as.character(counter), "_QC"))
      setwd(paste0(path_data_in, "/batch_", as.character(counter), "_QC"))
      
      if (plots == TRUE) {
        for (k in 1:dim(targExtracRes_QC[1])) {
          componentID = k
          
          plotdata_list_1 <- list()
          plotdata_list_std_1 <- list()
          for (s in 1:length(msRawData_QCs)) {
            componentuniqueID <-
              mergedata_list_QC[[s]][mergeData_QC$ID == k, ]$uniqueID
            componentname <- targExtracRes_QC$NAME[componentID]
            componentmz <- targExtracRes_QC$m.z[componentID]
            
            rawData <- msRawData_QCs[[s]]
            
            
            plotEICdata1 <-
              plotEIC(rawData, mzrange = mzRanges[componentuniqueID, ], rtrange = rtRanges[componentuniqueID, ])
            # plotdata1 <- as.data.frame(plotEICdata1@eic$xcmsRaw[[1]])
            colnames(plotEICdata1) <- cbind("rt", "intensity")
            plotdata1 <- as.data.frame(plotEICdata1)
            plotdata1$rt <- as.numeric(plotdata1$rt)
            plotdata1$rt <- plotdata1$rt / 60
            
            
            plotdata_list_1 <- c(plotdata_list_1, list(plotdata1))
            
            
          }
          
          
          p1 <-
            ggplot(
              data = bind_rows(plotdata_list_1, .id = "Sample"),
              aes(
                x = rt,
                y = intensity,
                colour = Sample
              )
            ) +
            scale_color_brewer(type = "qual", palette = "Set1") +
            geom_line() +
            geom_point() +
            ylab("Intensity") +
            xlab("RT (min)") +
            theme_light() +
            geom_vline(
              xintercept = targExtracRes_QC$tr[componentID] / 60,
              linetype = "dotted",
              size = 1,
              col = "red"
            ) +
            geom_vline(
              xintercept = targExtracRes_QC$trOfPeak[componentID] / 60,
              linetype = "dotted",
              size = 1,
              col = "blue"
            )
          
          
          rtminrange = rtRanges[componentuniqueID, 1] - 120
          rtmaxrange = rtRanges[componentuniqueID, 2] + 120
          if (rtminrange < min(rawData@scantime)) {
            rtminrange = min(rawData@scantime)
          }
          if (rtmaxrange > max(rawData@scantime)) {
            rtmaxrange = max(rawData@scantime)
          }
          
          plotdata_list_2 <- list()
          plotdata_list_std_2 <- list()
          for (s in 1:length(msRawData_QCs)) {
            rawData <- msRawData_QCs[[s]]
            plotEICdata2 <-
              plotEIC(
                rawData,
                mzrange = mzRanges[componentuniqueID, ],
                rtrange = cbind(rtminrange, rtmaxrange)
              )
            colnames(plotEICdata2) <- cbind("rt", "intensity")
            plotdata2 <- as.data.frame(plotEICdata2)
            plotdata2$rt <- as.numeric(plotdata2$rt)
            plotdata2$rt <- plotdata2$rt / 60
            
            
            plotdata_list_2 <- c(plotdata_list_2, list(plotdata2))
            
            
            
          }
          
          p2 <-
            ggplot(
              data = bind_rows(plotdata_list_2, .id = "Sample"),
              aes(
                x = rt,
                y = intensity,
                colour = Sample
              )
            ) +
            geom_line() +
            ylab("Intensity") +
            scale_color_brewer(type = "qual", palette = "Set1") +
            xlab('RT (min)') +
            scale_x_continuous(breaks = round(seq(
              min(plotdata2$rt), max(plotdata2$rt), by = 0.5
            ), 1)) +
            theme_light() +
            geom_vline(
              xintercept = targExtracRes_QC$tr[componentID] / 60,
              linetype = "dotted",
              size = 1,
              col = "red"
            )
          
          
          
          
          
          #plots samen plotten met passende titel etc.
          
          p <- p1 / p2
          p <- p + plot_annotation(
            title = paste(
              "Extracted EIC of",
              componentname,
              "/n m/z =",
              as.character(mzRanges[componentuniqueID, 1]),
              "-",
              as.character(mzRanges[componentuniqueID, 2]),
              "/n RT_db =",
              as.character(targExtracRes_QC$tr[componentID] / 60),
              "-",
              "Int=",
              as.character(targExtracRes_QC$peakHeight[componentID])
            )
            ,
            theme = theme_light()
          )
          file = paste0(componentname, ".png")
          file <- gsub(":", "", file)
          
          ggsave(file, plot = p)
          
        }
      }
      
      # Save intensity table of QC's
      write.csv(int_table_QC, "intensity_table_QC.csv")
      saveRDS(targExtracRes_list_QC, "targextraclist_QC.RDS")
      
      #If QC's ok --> change search RT to mean of RT where peaks where found in QC's
      
      #replace "search RT" with average of retention time found in the QC's
      retention_times <-
        simplify2array(targExtracRes_list_QC)["trOfPeak", ]
      retention_times_mean <-
        Reduce("+", retention_times) / length(retention_times)
      
      dbData$trold <- dbData$tr
      dbData$tr <- retention_times_mean
      
      
      
      
      #
      #if peak is NF, change RT back to original database RT to avoid errors
      for (i in 1:length(dbData$tr)) {
        if (is.na(dbData$tr[i]) == TRUE) {
          if (polarity == "positive") {
            dbData$tr[i] <- masslist_positive$tr[i]
          }
          else{
            dbData$tr[i] <- masslist_negative$tr[i]
          }
        }
      }
      
      
      #Re-run processing with updated RT
      
      start <- Sys.time()
      
      #insert raw data samples
      setwd(path_data_in)
      samplefiles <- batch_list[[counter]]
      xcmsset2 <- xcmsSet(msRawData_samples[samplefiles])
      msRawData_Samples <-
        getXcmsRaw(xcmsset2, sampleidx = seq(1, length(msRawData_samples[samplefiles]), 1))
      
      
      mergedata_list <- list()
      targExtracRes_list <- list()
      
      for (p in 1:length(msRawData_samples[samplefiles])) {
        #
        #Load raw data
        rawData <- msRawData_Samples[[p]]
        
        #select mz & rt from database
        mzmed <- dbData$`m/z`
        rtmed <- dbData$tr
        #combine m/z & rt & give unique identifier
        mzAndRt <-
          as.data.frame(cbind(c(1:length(mzmed)), paste(mzmed, rtmed)))
        colnames(mzAndRt) <- c("ID", "mzAndRt")
        mzAndRt$ID <- as.numeric(mzAndRt$ID)
        
        #find unique combinations of mz & RT and give identifier
        uniqueMzAndRtID <- which(!duplicated(mzAndRt$mzAndRt))
        uniqueMzAndRt <- mzAndRt[uniqueMzAndRtID, ]
        uniqueMzAndRt$ID <- c(1:nrow(uniqueMzAndRt))
        colnames(uniqueMzAndRt) <- c("uniqueID", "mzAndRt")
        
        #merge mz & RT combinations with their unique ID
        mergeData <- merge(mzAndRt, uniqueMzAndRt, by = "mzAndRt")
        mergeData <- mergeData[order(mergeData$ID), ]
        
        
        #save mz & rt from the unique mz & RT combo's
        mzmed <- mzmed[uniqueMzAndRtID]
        rtmed <- rtmed[uniqueMzAndRtID]
        
        
        #calculate delta m/z for each mass based on input ppm
        mzdeltas <- sapply(mzmed, function(mzmed)
          mzmed * ppm / 10 ^ 6)
        
        #calculate mzrange based on delta mz
        mzRanges <-
          cbind(as.numeric(mzmed) - mzdeltas,
                as.numeric(mzmed) + mzdeltas)
        
        #if an upper m/z boundary is lower than the minimum m/z range, set it to the minimum m/z
        indexTemp <- which(mzRanges[, 2] < min(rawData@mzrange))
        mzRanges[indexTemp, ] <- min(rawData@mzrange)
        
        #if a lower m/z boundary is higher than the max m/z, set it to the max m/z
        indexTemp <- which(mzRanges[, 1] > max(rawData@mzrange))
        mzRanges[indexTemp, ] <- max(rawData@mzrange)
        
        #if an upper limit is higher than the max mz & the lower limit is smaller than the max m/z, set the upper limit to the max m/z
        mzRanges[which(mzRanges[, 2] > max(rawData@mzrange) &
                         mzRanges[, 1] < max(rawData@mzrange)), 2] <- max(rawData@mzrange)
        
        #if a upper limit is larger than the minimum & the lower limit is lower than the minimum, set the lower limit to the min m/z
        mzRanges[which(mzRanges[, 2] > min(rawData@mzrange) &
                         mzRanges[, 1] < min(rawData@mzrange)), 1] <- min(rawData@mzrange)
        
        
        #calculate rt range based on deltaTR
        rtRanges <-
          cbind(as.numeric(rtmed) - deltaTR / 2,
                as.numeric(rtmed) + deltaTR / 2)
        
        #for al upper RT limits -2 lower than minimum RT, change to lower limit of min rt & upper limit to min + 10
        indexTemp <- which(rtRanges[, 2] - 2 < min(rawData@scantime))
        rtRanges[indexTemp, ] <-
          cbind(rep(min(rawData@scantime), length(indexTemp)), rep(min(rawData@scantime) +
                                                                     10, length(indexTemp)))
        
        
        #for al lower RT limits +2 higher than maximum RT, change lower limit to max -10 and upper limit to max
        indexTemp <- which(rtRanges[, 1] + 2 > max(rawData@scantime))
        rtRanges[indexTemp, ] <-
          cbind(rep(max(rawData@scantime) - 10, length(indexTemp)), rep(max(rawData@scantime), length(indexTemp)))
        
        
        #for upper limits higher than max and lower limit lower than max, change upper limit to max
        rtRanges[which(rtRanges[, 2] > max(rawData@scantime) &
                         rtRanges[, 1] < max(rawData@scantime)), 2] <- max(rawData@scantime)
        #for upper limits higher than min and lower limit lower than min, change lower limit to min
        rtRanges[which(rtRanges[, 2] > min(rawData@scantime) &
                         rtRanges[, 1] < min(rawData@scantime)), 1] <- min(rawData@scantime)
        
        #extract EIC data for defined ranges
        # EICdata <- getEIC(rawData, mzrange = mzRanges, rtrange = rtRanges) #this function uses profile matrix data
        
        rawEICdata_list <- list()
        for (i in 1:dim(rtRanges)[1]) {
          rawEICdata <-
            plotEIC(rawData, mzrange = mzRanges[i, ], rtrange = rtRanges[i, ]) #using this we can extract raw data
          rawEICdata_list[[i]] <- rawEICdata
        }
        
        
        #extract raw data
        # EICdataEIC <- EICdata@eic$xcmsRaw
        # #select the ones that have a uniqueID
        # EICdataEIC <- EICdataEIC[mergeData$uniqueID]
        rawEICdata_list <- rawEICdata_list[mergeData$uniqueID]
        
        #list to char function to create vectors with all EIC data per mzrange & rtrange
        list2Character <- function(ithMatrix) {
          return(paste(paste(ithMatrix[, 1], ithMatrix[, 2], sep = " "), collapse = ";"))
        }
        
        
        
        # EICdataEICVector <- unlist(lapply(EICdataEIC, list2Character))
        rawEICdatavector <-
          unlist(lapply(rawEICdata_list, list2Character))
        #add to dbdata
        dbData$EIC <- rawEICdatavector
        
        #functie maken die "peak detection" en entrophy calc toepast op elke set EIC.
        # binnen elke EIC lokale maxima selecteren en als output RT & hoogte van lokaal maximum
        #er zit geen functie voor area te berekenen in de functie
        #berekend entropie binnen trrange interval
        func <- function(ithRowDbData, ithEICdataEIC) {
          extractedPeaks <-
            peakDectAndEntroCal(ithEICdataEIC, trRange = trRange, m = m)
          return(cbind(ithRowDbData[rep(1, nrow(extractedPeaks)), ], extractedPeaks))
        }
        
        dbDataList <- split(dbData, 1:nrow(dbData))
        
        targExtracRes <-
          mapply(func, dbDataList, rawEICdata_list, SIMPLIFY = F)
        targExtracRes <-
          do.call(rbind, lapply(targExtracRes, data.frame))
        targExtracRes$trOfPeak <- as.numeric(targExtracRes$trOfPeak)
        
        
        
        mergedata_list <- c(mergedata_list, list(mergeData))
        targExtracRes_list <- c(targExtracRes_list, list(targExtracRes))
        
      }
      
      int_table <- data.frame(targExtracRes_list[[1]]$NAME)
      
      for (k in 1:length(targExtracRes_list)) {
        int_table <- cbind(int_table, targExtracRes_list[[k]]$peakHeight)
      }
      
      colnames(int_table) <-
        append("Name", msRawData_samples[samplefiles])
      
      
      
      
      dir.create(paste0(
        path_data_in,
        "/batch_",
        as.character(counter),
        "_samples"
      ))
      setwd(paste0(
        path_data_in,
        "/batch_",
        as.character(counter),
        "_samples"
      ))
      
      
      if(plots == TRUE){
      #data plotten met "extracted" range en aangepaste range --> hier plots met "getEIC"
      for (k in 1:dim(targExtracRes[1])) {
        componentID = k
        
        plotdata_list_1 <- list()
        
        for (s in 1:length(msRawData_samples[samplefiles])) {
          componentuniqueID <-
            mergedata_list[[s]][mergeData$ID == k, ]$uniqueID
          componentname <- targExtracRes$NAME[componentID]
          componentmz <- targExtracRes$m.z[componentID]
          
          rawData <- msRawData_Samples[[s]]
          
          
          plotEICdata1 <-
            plotEIC(rawData, mzrange = mzRanges[componentuniqueID, ], rtrange = rtRanges[componentuniqueID, ])
          # plotdata1 <- as.data.frame(plotEICdata1@eic$xcmsRaw[[1]])
          colnames(plotEICdata1) <- cbind("rt", "intensity")
          plotdata1 <- as.data.frame(plotEICdata1)
          plotdata1$rt <- as.numeric(plotdata1$rt)
          plotdata1$rt <- plotdata1$rt / 60
          
          
          plotdata_list_1 <- c(plotdata_list_1, list(plotdata1))
          
          
        }
        
        
        p1 <-
          ggplot(data = bind_rows(plotdata_list_1, .id = "Sample"),
                 aes(
                   x = rt,
                   y = intensity,
                   colour = Sample
                 )) +
          scale_color_grey() +
          geom_line() +
          geom_point() +
          ylab("Intensity") +
          xlab("RT (min)") +
          theme_light() +
          theme(legend.position = "none") +
          geom_vline(
            xintercept = targExtracRes$tr[componentID] / 60,
            linetype = "dotted",
            size = 1,
            col = "red"
          )
        #geom_vline(xintercept = targExtracRes$trOfPeak[componentID]/60, linetype = "dotted", size = 1, col = "blue")
        
        
        rtminrange = rtRanges[componentuniqueID, 1] - 120
        rtmaxrange = rtRanges[componentuniqueID, 2] + 120
        if (rtminrange < min(rawData@scantime)) {
          rtminrange = min(rawData@scantime)
        }
        if (rtmaxrange > max(rawData@scantime)) {
          rtmaxrange = max(rawData@scantime)
        }
        
        plotdata_list_2 <- list()
        plotdata_list_std_2 <- list()
        for (s in 1:length(msRawData_samples[samplefiles])) {
          rawData <- msRawData_Samples[[s]]
          plotEICdata2 <-
            plotEIC(rawData,
                    mzrange = mzRanges[componentuniqueID, ],
                    rtrange = cbind(rtminrange, rtmaxrange))
          # plotdata1 <- as.data.frame(plotEICdata1@eic$xcmsRaw[[1]])
          colnames(plotEICdata2) <- cbind("rt", "intensity")
          plotdata2 <- as.data.frame(plotEICdata2)
          plotdata2$rt <- as.numeric(plotdata2$rt)
          plotdata2$rt <- plotdata2$rt / 60
          
          
          plotdata_list_2 <- c(plotdata_list_2, list(plotdata2))
          
          
          
        }
        
        p2 <-
          ggplot(data = bind_rows(plotdata_list_2, .id = "Sample"),
                 aes(
                   x = rt,
                   y = intensity,
                   colour = Sample
                 )) +
          geom_line() +
          ylab("Intensity") +
          scale_color_grey() +
          xlab('RT (min)') +
          #geom_line(data = bind_rows(plotdata_list_std_2, .id = "Sample"),aes(x=rt,y=intensity,colour = Sample),linetype = "dashed",size=1) +
          scale_x_continuous(breaks = round(seq(
            min(plotdata2$rt), max(plotdata2$rt), by = 0.5
          ), 1)) +
          theme_light() +
          theme(legend.position = "none") +
          geom_vline(
            xintercept = targExtracRes$tr[componentID] / 60,
            linetype = "dotted",
            size = 1,
            col = "red"
          )
        #geom_vline(xintercept = targExtracRes$trOfPeak[componentID]/60, linetype = "dotted", size = 1, col = "blue")
        
        
        
        
        #plots samen plotten met passende titel etc.
        
        p <- p1 / p2
        p <- p + plot_annotation(
          title = paste(
            "Extracted EIC of",
            componentname,
            "/n m/z =",
            as.character(mzRanges[componentuniqueID, 1]),
            "-",
            as.character(mzRanges[componentuniqueID, 2]),
            "/n RT_db =",
            as.character(targExtracRes$tr[componentID] / 60),
            "-",
            "Int=",
            as.character(targExtracRes$peakHeight[componentID])
          )
          ,
          theme = theme_light()
        )
        file = paste0(componentname, ".png")
        file <- gsub(":", "", file)
        
        ggsave(file, plot = p)
        
      }
      }
      
      write.csv(int_table, "intensity_table.csv")
      saveRDS(targExtracRes_list, "targextraclist.RDS")
      stop <- Sys.time()
      runtime <- stop - start
      print(runtime)
      
    }
  }