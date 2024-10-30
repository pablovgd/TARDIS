#' Sum Intensities of Spectra


.sum_intensities <- function(x, ...) {
  if (nrow(x)) {
    cbind(mz = NA_real_,
          intensity = sum(x[, "intensity"], na.rm = TRUE))
  } else {
    cbind(mz = NA_real_, intensity = NA_real_)
  }
}
