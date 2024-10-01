#' Check if any samples are missing spectra
#'
#' @param spectra `Spectra` object
#'
#' @importFrom ProtGenerics dataOrigin
#' @importFrom Spectra isEmpty
#' @return
#' @export
#'
#' @examples
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
