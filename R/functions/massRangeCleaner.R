#' Cleaning mass ranges for processing spectral data
#'
#' As several functions within the peak detection framework do not cope well with, this function selects a mass window from the files and writes new files back to the input path.
#'
#' @param path Path to .mzML or .mzXML that contain multiple mass window scan ranges
#' @param massrange Mass range of interest. Can only contain masses from one mass window
#'
#' @return
#' @export
#'
#' @examples
massRangeCleaner <- function(path, massrange){
  
  
  datafiles <- list.files(path = path, pattern = ".mzML|.mzXML", recursive = TRUE)
  lipid_data <- readMSData(datafiles,msLevel. = 1,mode = "onDisk")
  lipid_data <- filterMz(lipid_data,mz= massrange)
  lipid_data <- filterEmptySpectra(lipid_data)
  writeMSData(lipid_data,paste("cleaned",datafiles,sep="_"))
  
}