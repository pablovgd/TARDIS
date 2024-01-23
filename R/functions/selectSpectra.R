selectSpectra <- function(path, QCPattern, SamplePattern){


datafiles <- list.files(path = path, pattern = ".mzML|.mzXML", recursive = TRUE)
msRawData_samples <-  datafiles[grep(pattern = SamplePattern,datafiles)]
msRawData_QC <-  datafiles[grep(pattern = QCPattern,datafiles)]


returnlist <- list(msRawData_samples,msRawData_QC)

return(returnlist)
}

