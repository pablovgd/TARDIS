#' @title  Check if any samples are missing spectra
#' @description
#' Display error if they do.
#' Currently stops if a sample is detected that has less than 90% of the mean
#' of spectra in all samples.
#'
#'
#' @param spectra `Spectra` object
#'
#' @importFrom ProtGenerics dataOrigin
#' @export
#' @author Pablo Vangeenderhuysen
checkScans <- function(spectra){
  scans_per_sample <- table(dataOrigin(spectra))
  mean <- (mean(scans_per_sample))
  bad_runs <- which(scans_per_sample < 0.9*mean)
  if(isEmpty(bad_runs) == FALSE){
    names <- basename(names(bad_runs))
    stop(paste("File",names,"contains less than 50% of the mean of scans
                   in the samples."))
  }
}
