selectSpectra <- function(path, QCPattern, SamplePattern){


datafiles <- list.files(path = path, pattern = ".mzML|.mzXML", recursive = TRUE)
msRawData_samples <-  datafiles[grep(pattern = SamplePattern,datafiles)]
msRawData_QC <-  datafiles[grep(pattern = ,datafiles)]




return(msRawData_samples,msRawData_QC)
}

