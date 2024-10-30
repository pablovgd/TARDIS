#' qscoreCalculator
#'
#' Implementation of the work by William Kumler to calculate quality metrics:
#' https://github.com/wkumler/MS_metrics
#'
#' @param rt
#' @param int


qscoreCalculator <- function(rt, int){
  #Check for bogus EICs
  if (length(rt)<5) {
    return(list(SNR = NA_real_, peak_cor = NA_real_))
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (rt - min(rt)) / (max(rt) - min(rt))

  # Create a couple different skews and test fit
  maybe_skews <- c(2.5, 3, 4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, int)


  #Calculate the normalized residuals
  residuals <- int / max(int) - perf_peak / max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while (new_res_sd < old_res_sd) {
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(int) - min(int)) / sd(norm_residuals * max(int))
  #Return the quality score
  return(list(SNR = SNR, peak_cor = peak_cor))
}
