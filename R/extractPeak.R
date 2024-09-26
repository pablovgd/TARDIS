#' filter Spectra to single peak in single sample
filterSpectra <- function(spectra,dataOrigin,rt_range,mz_range) {
  spectra <- spectra |>
    filterDataOrigin(dataOrigin) |>
    filterRt(rt_range) |>
    filterMzRange(mz_range)
  spectra
}

#' Sum Intensities of Spectra
.sum_intensities <- function(x, ...) {
  if (nrow(x)) {
    cbind(mz = NA_real_,
          intensity = sum(x[, "intensity"], na.rm = TRUE))
  } else
    cbind(mz = NA_real_, intensity = NA_real_)
}

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

