#' Check if any samples are missing more than 10% of scans (compared to the
#' average number of scans across samples)
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
  mean_scans <- (mean(scans_per_sample))
  bad_runs <- which(scans_per_sample < (0.9 * mean_scans))
  if (isEmpty(bad_runs) == FALSE) {
    names <- basename(names(bad_runs))
    stop(paste("File", names,
               "contains less than 90% of the mean of scans across samples.\n"))
  }
}
