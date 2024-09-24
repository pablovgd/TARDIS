#' @title Function to find peaks given target info and spectral data
#' 
#' @param spectra a `Spectra` object.
#' @param smoothing `logical` `TRUE` or `FALSE`.
#' 
#' @importFrom signal sgolayfilt
#' 
#' 
#' @author Pablo Vangeenderhuysen
#' 
#' @return named `numeric(4)` with x: left and right retention time borders
#' y: left and right intensity borders, rt: retention vector of peak, int:
#' intensity vector of peak.
#'    
#' @export
#' 
find_peak <- function(spectra,smoothing,search_rt){
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
  border <- find_peak_points(rt, smoothed, search_rt,
                             .check = FALSE)
  border
}