#' @title  Check if any samples are missing spectra
#' @description
#' Display warning if they do.
#' Currently warns user if a sample is detected that has less than 50% of the mean
#' of spectra in all samples.
#'
#' @param spectra `Spectra` object
#'
#' @importFrom ProtGenerics dataOrigin
#' @importFrom Spectra isEmpty
#' @export
#' @author Pablo Vangeenderhuysen
checkScans <- function(spectra){
  scans_per_sample <- table(dataOrigin(spectra))
  mean <- (mean(scans_per_sample))
  bad_runs <- which(scans_per_sample < 0.5*mean)
  if(isEmpty(bad_runs) == FALSE){
    names <- basename(names(bad_runs))
    warning(paste("File",names,"contains less than 50% of the mean of scans
                   in the samples."))
  }
}
