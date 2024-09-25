#' @title Function to extract EIC from Spectra object
#' 
#' @param spectra a `Spectra` object.
#' 
#' 
#' 
#' @author Pablo Vangeenderhuysen
#' 
#' @return
#'    
#' @export
#' 
extract_eic <- function(spectra){
  sfs_agg <-
    addProcessing(sample_spectra, .sum_intensities)
  eic <-
    cbind(rtime(sfs_agg),
          unlist(intensity(sfs_agg), use.names = FALSE))
  rownames(eic) <- NULL
  colnames(eic) <- c("rt","int")
  eic
}
