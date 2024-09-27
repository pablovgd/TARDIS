#' @title T.A.R.D.I.S. Peak Detection
#'
#' @description
#' Main function of the T.A.R.D.I.S. package that is called in the Shiny app.
#' Given data files and a list of targeted compounds it returns the area of
#' those peaks, optional diagnostic plots and several other parameters.
#' See vignette for a detailed tutorial.
#'
#' @param file_path `character(1)`: path to the .mzML or .mzXML files
#'     containing LC-MS data.
#' @param dbData Output of [createTargetList()]
#' @param ppm `numeric(1)` Allowed deviance from given m/z of targets in ppm.
#' @param rtdev `numeric(1)` Allowed deviance from given retention time of
#'     compound, defines search window for the peak picking algorithm.
#' @param mass_range `numeric(2)` If the user uses data with overlapping mass
#'     windows, only one mass window at the time can be analyzed, specify this
#'     here.
#' @param polarity `character(1)` Ionisation mode to be considered, can be
#'     either "positive" or "negative"
#' @param output_directory `character(1)` Provide directory to store output
#' @param plots_samples `logical(1)` Create plots for all samples
#' @param plots_QC `logical(1)` Create plots for all QCs
#' @param diagnostic_plots `logical(1)` Create diagnostic plots of 5 QCs
#'     spread across the runs
#' @param batch_positions `list` Indicate start and end file of each batch,
#'     e.g. `list(c(1,20),c(21,40))`
#' @param QC_pattern `character(1)`  Pattern of QC files
#' @param sample_pattern `character(1)` Pattern of sample files
#' @param rt_alignment `logical(1)` Align retention time based on internal
#'     standard compounds in the QC samples.
#' @param int_std_id `character` Provide ID's of internal standard compounds for
#'     retention time alignment
#' @param screening_mode `logical(1)` Run the algorithm over 5 QCs to quickly
#'     check retention time shifts
#' @param smoothing `logical(1)` Smooth the peaks with [sgolayfilt()],
#'
#' @import MsExperiment
#' @importFrom Spectra MsBackendMzR
#' @importFrom Spectra filterMzRange
#' @importFrom Spectra filterEmptySpectra
#' @importFrom Spectra filterDataOrigin
#' @importFrom Spectra filterRt
#' @importFrom Spectra dataOrigin
#' @importFrom Spectra addProcessing
#' @importFrom signal sgolayfilt
#' @importFrom xcms PeakGroupsParam
#' @importFrom xcms adjustRtime
#' @importFrom xcms applyAdjustedRtime
#' @importFrom xcms rtime
#' @importFrom xcms intensity
#' @importFrom pracma trapz
#' @importFrom BiocParallel SnowParam
#' @importFrom tidyr spread
#' @importFrom writexl write_xlsx
#' @importFrom dplyr summarise
#' @importFrom dplyr summarise_at
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr first
#'
#' @return
#' @export
#'

## jo: wouldn't it be better to call the function on a data object instead
## of a file path? The (advanced) user could eventually do some more quality
## checks on the data before?
## Pablo: Definitely something I should do, but I might do it a separate
## function, since for use with the GUI the file input is nice.

tardis_peaks <-
  function(file_path,
           dbData,
           ppm = 5,
           rtdev = 18,
           mass_range = c(67,1000),
           polarity = "positive",
           output_directory,
           plots_samples = FALSE,
           plots_QC = FALSE,
           diagnostic_plots = TRUE,
           batch_positions,
           QC_pattern = "QC",
           sample_pattern = "",
           rt_alignment = TRUE,
           int_std_id,
           screening_mode = FALSE,
           smoothing = TRUE) {
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
    files <-
      list.files(file_path, full.names = T, pattern = "mzML|mzXML")
    if (is.null(mass_range) == FALSE) {
      dbData <-
        dbData[which(dbData$`m/z` < mass_range[2] &
                       dbData$`m/z` > mass_range[1]), ]
    }
    info_compounds <- dbData
    if (screening_mode == TRUE) {
      results_screening <-
        data.frame(
          Component = character(0),
          Sample = character(0),
          AUC = numeric(0),
          SNR = numeric(0),
          peak_cor = numeric(0),
          foundRT = numeric(0),
          pop = numeric(0)
        )
      QC_files <-
        files[grep(pattern = QC_pattern, files)]
      data_QC <- MsExperiment::readMsExperiment(
        spectraFiles = QC_files,
        backend = MsBackendMzR(),
        BPPARAM = SnowParam(workers = 1L)
      )
      if (is.null(mass_range) == FALSE) {
        spectra_QC <- data_QC@spectra |>
          filterMzRange(mass_range) |>
          filterEmptySpectra()
      } else{
        spectra_QC <- data_QC@spectra
        }
        checkScans(spectra_QC)

      ## Create ranges for all compounds
      ranges <- createRanges(data_QC, dbData, ppm, rtdev)
      ## Get mz & rt ranges
      mzRanges <- ranges[[1L]]
      rtRanges <- ranges[[2L]]

      if (rt_alignment == TRUE) {
        ## Get the ranges for the internal standard compounds
        internal_standards_rt <-
          rtRanges[which(dbData$ID %in% int_std_id), ]
        internal_standards_mz <-
          mzRanges[which(dbData$ID %in% int_std_id), ]
        dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
        ## Get QC sample names
        sample_names <-
          lapply(data_QC@sampleData$spectraOrigin, basename)
        ## Initiate empty vectors
        int_std_foundrt <- c()
        int_std <- c()
        ## Retrieve foundRT of internal standards in QC's,
        ## loop over all samples and all internal standards
        for (j in 1:dim(internal_standards_rt)[1]) {
          rt_list = list()
          int_list = list()
          x_list = list()
          y_list = list()
          for (i in 1:length(sample_names)) {
            sample_name <- unlist(sample_names[i])
            filtered_spectra <- filterSpectra(spectra_QC,
                                             unique(dataOrigin(spectra_QC))[i],
                                             internal_standards_rt[j, ],
                                             internal_standards_mz[j, ])
            eic <- extract_eic(filtered_spectra)
            rt <- eic[, 1L]
            int <- eic[, 2L]
            #NA intensities are set to zero --> should change this so only NA's
            #at the edges get changed to zero, so the ones IN the peak will be
            #imputed
            int[which(is.na(int))] = 0
            #To determine the borders more easily, smoothing is applied
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            if (smoothing == TRUE) {
              int <- smoothed
              int[int < 0] <- 0
            }
            #Border detection
            border <-
              find_peak_points(rt, smoothed, dbData_std$tr[j],
                               .check = FALSE)

            #Save found RT for internal standard target
            int_std_foundrt <-
              cbind(int_std_foundrt, rt[border[3L]]) #this will finally contain
              #all found rts for the different internal standards in this sample
          }

          ## this will contain all the found rts in all samples
          int_std <-
            rbind(int_std, int_std_foundrt)
          int_std_foundrt <- c()
        }
        ## Define parameters for retention time adjustment, based on the QC's
        ## and the internal standards
        param <-
          PeakGroupsParam(
            minFraction = 0.9,
            span = 0.5  ,
            peakGroupsMatrix = int_std
          )
        data_QC <- adjustRtime(data_QC, param = param)
        data_QC <- applyAdjustedRtime(data_QC)
      }
      #Find all targets in x QC's
      if (is.null(mass_range) == FALSE) {
        spectra_QC <- data_QC@spectra |>
         filterMzRange(mass_range)
         filterEmptySpectra()
      } else{
        spectra_QC <- data_QC@spectra
      }
      for (j in 1:dim(rtRanges)[1]) {
        rt_list = list()
        int_list = list()
        x_list = list()
        y_list = list()
        for (i in 1:length(sample_names)) {
          sample_name <- unlist(sample_names[i])
          filtered_spectra <- filterSpectra(spectra_QC,
                                            unique(dataOrigin(spectra_QC))[i],
                                            rtRanges[j, ],
                                            mzRanges[j, ])
          eic <- extract_eic(filtered_spectra)
          rt <- eic[, 1L]
          int <- eic[, 2L]
          int[which(is.na(int))] = 0
          smoothed <- sgolayfilt(int, p = 3, n = 7)
          if (smoothing == TRUE) {
            int <- smoothed
            int[int < 0] <- 0
          }
          border <-
            find_peak_points(rt, smoothed, dbData$tr[j], .check = FALSE)
          idx <- border[1L]:border[2L]
          x <- rt[idx]
          y <- int[idx]
          rt_list <- c(rt_list, list(rt))
          int_list <- c(int_list, list(int))
          x_list <- c(x_list, list(x))
          y_list <- c(y_list, list(y))
          ## Check if there are at least two unique values for the component
          if (length(unique(y)) > 1) {
            auc <- trapz(x, y)
            pop <- length(x)
            qscore <- qscoreCalculator(x, y)
            compound_info <- dbData[j, ]
            results_screening <- rbind(
              results_screening,
              data.frame(
                Component = compound_info$ID,
                Sample = sample_name,
                AUC = auc,
                MaxInt = int[border[3L]],
                SNR = qscore[1],
                peak_cor = qscore[2],
                foundRT = rt[border[3L]],
                pop = pop,
                compound_info
              )
            )
          } else{
            compound_info <- dbData[j, ]
            # Append results to the data frame with compound information
            results_screening <- rbind(
              results_screening,
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
      avg_metrics_table <- results_screening %>%
        group_by(Component) %>%
        summarise_at(vars(-Sample), list(~ if (is.numeric(.))
          mean(., na.rm = TRUE)
          else
            first(.)))
      write.csv(avg_metrics_table,
                file = paste0(output_directory, "qc_screening.csv"))
    } else {
      ## Loop over the batches
      for (batchnr in 1:length(batch_positions)) {
        dbData <- info_compounds #need to reset? better to keep updated from
        # last batch?
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
        files_batch <-
          files[batch_positions[[batchnr]][1]:batch_positions[[batchnr]][2]]
        data_batch <- readMsExperiment(
          spectraFiles = files_batch,
          backend = MsBackendMzR(),
          BPPARAM = SnowParam(workers = 1)
        )
        checkScans(data_batch@spectra)
        #Define study and QC samples --> all not QC files are deemed study files
        sampleData(data_batch)$sample_type <- "study"
        sampleData(data_batch)$sample_type[grep(pattern = QC_pattern,
                                                files_batch)] <- "QC"
        data_QC <-
          data_batch[which(sampleData(data_batch)$sample_type == "QC")]
        if (is.null(mass_range == FALSE)) {
          spectra_QC <- data_QC@spectra |>
            filterMzRange(mass_range) |>
            filterEmptySpectra()
        } else{
          spectra_QC <- data_QC@spectra
        }
        ranges <- createRanges(data_QC, dbData, ppm, rtdev)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        if (rt_alignment == TRUE) {
          ## Get the ranges for the internal standard compounds
          internal_standards_rt <-
            rtRanges[which(dbData$ID %in% int_std_id), ]
          internal_standards_mz <-
            mzRanges[which(dbData$ID %in% int_std_id), ]
          dbData_std <- dbData[which(dbData$ID %in% int_std_id), ]
          sample_names <-
            lapply(data_QC@sampleData$spectraOrigin, basename)
          int_std_foundrt <- c()
          int_std <- c()
          for (j in 1:dim(internal_standards_rt)[1]) {
            rt_list = list()
            int_list = list()
            x_list = list()
            y_list = list()
            for (i in 1:length(sample_names)) {
              sample_name <- unlist(sample_names[i])
              spectra_filtered <- filterSpectra(spectra_QC,
                unique(dataOrigin(spectra_QC))[i],
                internal_standards_rt[j, ],
                internal_standards_mz[j, ])
              eic <- extract_eic(spectra_filtered)
              rt <- eic[, 1L]
              int <- eic[, 2L]
              int[which(is.na(int))] = 0
              smoothed <- sgolayfilt(int, p = 3, n = 7)
              if (smoothing == TRUE) {
                int <- smoothed
                int[int < 0] <- 0
              }
              border <-
                find_peak_points(rt, smoothed, dbData_std$tr[j],
                                 .check = FALSE)
              int_std_foundrt <-
                cbind(int_std_foundrt, rt[border[3L]])
            }
            int_std <-
              rbind(int_std, int_std_foundrt)
            int_std_foundrt <- c()
          }
          param <-
            PeakGroupsParam(
              minFraction = 0.9,
              span = 0.5  ,
              peakGroupsMatrix = int_std,
              subset = which(sampleData(data_batch)$sample_type == "QC")
            )
          data_batch <- adjustRtime(data_batch, param = param)
          data_batch <- applyAdjustedRtime(data_batch)
        }
        ## Now, we try and find ALL compounds in the QC samples and save their
        ## foundRT to search the compounds at that RT in the sample files
        ## Skip this step if there aren't any QC's available.
        data_QC <-
          data_batch[which(sampleData(data_batch)$sample_type == "QC")]
        if (length(data_QC) != 0) {
          sample_names <-
            lapply(data_QC@sampleData$spectraOrigin, basename)
          if (is.null(mass_range) == FALSE) {
            spectra_QC <- data_QC@spectra |>
              filterMzRange(mass_range) |>
              filterEmptySpectra()
          } else{
            spectra_QC <- data_QC@spectra
          }
          for (j in 1:dim(rtRanges)[1]) {
            rt_list = list()
            int_list = list()
            x_list = list()
            y_list = list()
            for (i in 1:length(sample_names)) {
              sample_name <- unlist(sample_names[i])
              filtered_spectra <- filterSpectra(spectra_QC,
                                          unique(dataOrigin(spectra_QC))[i],
                                          rtRanges[j, ],
                                          mzRanges[j,])
              eic <- extract_eic(filtered_spectra)
              rt <- eic[, 1L]
              int <- eic[, 2L]
              int[which(is.na(int))] = 0
              smoothed <- sgolayfilt(int, p = 3, n = 7)
              if (smoothing == TRUE) {
                int <- smoothed
                int[int < 0] <- 0
              }
              border <- find_peak_points(rt, smoothed, dbData$tr[j],
                                         .check = FALSE)
              idx <- border[1L]:border[2L]
              x <- rt[idx]
              y <- int[idx]
              rt_list <- c(rt_list, list(rt))
              int_list <- c(int_list, list(int))
              x_list <- c(x_list, list(x))
              y_list <- c(y_list, list(y))
              if (length(unique(y)) > 1) {
                auc <- trapz(x, y)
                pop <- length(x)
                qscore <- qscoreCalculator(x, y)
                compound_info <- dbData[j, ]
                results_QCs_batch <- rbind(
                  results_QCs_batch,
                  data.frame(
                    Component = compound_info$ID,
                    Sample = sample_name,
                    AUC = auc,
                    MaxInt = int[border[3L]],
                    SNR = qscore[1],
                    peak_cor = qscore[2],
                    foundRT = rt[border[3L]],
                    pop = pop,
                    compound_info
                  )
                )

              } else{
                compound_info <- dbData[j, ]
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
          ## Replace rtmed with average foundRT from previous results
          new_rt_avg <- results_QCs_batch %>%
            group_by(ID) %>%
            summarise(mean = mean(foundRT), na.rm = TRUE)
          dbData <- merge(dbData, new_rt_avg, by = "ID")
          dbData$trold <- dbData$tr
          dbData$tr <- new_rt_avg$mean
          ## If no RT is found, restore old RT
          for (k in 1:dim(dbData[1])) {
            if (is.na(dbData$tr[k]) == TRUE) {
              dbData$tr[k] <- dbData$trold[k]
            }
          }
        }
        ## Next do the whole analysis for the samples in the same batch of the
        ## QC's to find ALL the compounds at the corrected RT. (SAMPLES + QC)
        ## Get sample data
        sample_names <-
          lapply(data_batch@sampleData$spectraOrigin, basename)
        #Create ranges around new RT
        ranges <- createRanges(data_batch, dbData, ppm, rtdev)
        mzRanges <- ranges[[1]]
        rtRanges <- ranges[[2]]
        if (is.null(mass_range) == FALSE) {
          spectra <- data_batch@spectra |>
            filterMzRange(mz = mass_range) |>
            filterEmptySpectra()
        } else{
          spectra <- data_batch@spectra
        }

        for (j in 1:dim(rtRanges)[1]) {
          rt_list = list()
          int_list = list()
          x_list = list()
          y_list = list()
          for (i in 1:length(sample_names)) {
            sample_name <- unlist(sample_names[i])
            spectra_filtered <- filterSpectra(spectra,
                                              unique(dataOrigin(spectra))[i],
                                              rtRanges[j, ],
                                              mzRanges[j, ])
            eic <- extract_eic(spectra_filtered)
            eic <- eic[which(duplicated(eic[, 1]) == FALSE), ]#why is this here?
            rt <- eic[, 1L]
            int <- eic[, 2L]
            int[which(is.na(int))] = 0
            smoothed <- sgolayfilt(int, p = 3, n = 7)
            if (smoothing == TRUE) {
              int <- smoothed
              int[int < 0] <- 0
            }
            border <- find_peak_points(rt, smoothed, dbData$tr[j],
                                       .check = FALSE)
            idx <- border[1L]:border[2L]
            x <- rt[idx]
            y <- int[idx]
            rt_list <- c(rt_list, list(rt))
            int_list <- c(int_list, list(int))
            x_list <- c(x_list, list(x))
            y_list <- c(y_list, list(y))
            # Check if there are at least two unique values for the component
            if (length(unique(y)) > 1) {
              auc <- trapz(x, y)
              pop <- length(x)
              qscore <- qscoreCalculator(x, y)
              compound_info <- dbData[j, ]
              results_samples <- rbind(
                results_samples,
                data.frame(
                  Component = compound_info$ID,
                  Sample = sample_name,
                  AUC = auc,
                  MaxInt = int[border[3L]],
                  SNR = qscore[1],
                  peak_cor = qscore[2],
                  foundRT = rt[border[3L]],
                  pop = pop,
                  compound_info
                )
              )

            } else{
              compound_info <- dbData[j, ]
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
      write.csv(peakcor_table,
                file = paste0(output_directory, "peakcor_table.csv"))
      #summarize feature table based on QC's
      avg_metrics_table <- NULL
      if (length(data_QC) != 0) {
        QC_results <-  results[grep("QC", results$Sample), ]
        avg_metrics_table <- QC_results %>%
          group_by(Component) %>%
          summarise_at(vars(-Sample), list(~ if (is.numeric(.))
            mean(., na.rm = TRUE)
            else
              first(.)))
        write_xlsx(avg_metrics_table,paste0(output_directory,
                                            "feat_table.xlsx"))
      }
      return(list(auc_table, avg_metrics_table))
    }
  }
