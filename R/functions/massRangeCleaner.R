massRangeCleaner <- function(path, massrange){
  
  
  datafiles <- list.files(path = path, pattern = ".mzML|.mzXML", recursive = TRUE)
  lipid_data <- readMSData(datafiles,msLevel. = 1,mode = "onDisk")
  lipid_data <- filterMz(lipid_data,mz= massrange)
  lipid_data <- filterEmptySpectra(lipid_data)
  writeMSData(lipid_data,paste("cleaned",datafiles,sep="_"))
  
}