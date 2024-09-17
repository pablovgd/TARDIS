#' Helper function to find the borders of a peak in a sequence of intensities
#'
#' @param sign_change `logical` with positions of local minima and maxima.
#'     Length of `sign_change` needs to be the same than the length of the
#'     original intensity vector. Eventually add leading and trailing `FALSE`.
#'
#' @param index `integer(1)` with the index of the detected peak
#'
#' @param min_dist `integer(1)` defining the minimum distance of the left and
#'     right border to the apex position. With `min_dist = 2` (the default)
#'     the right border needs to be >= `index + 2`
#'
#' @return named `numeric(2)` with the left and right index of the border. If
#'     no local minima were found 1 and `length(sign_change)` are returned.
#' 
#' @author Pablo Vangeenderhuysen
#' 
#' @noRd
.find_peak_border <- function(sign_change, index, min_dist = 2) {
    l <- length(sign_change)
    left_true <- 1L
    right_true <- l
    if (index == 1)
        index <- index + 1
    ## Search for TRUE to the left --> if sign changes and at least doens't
    ## change at the next point
    for (i in (index - 1):1) {
        if (sign_change[i] && (index - i) >= min_dist) {
            left_true <- i
            break
        }
    }
    if (index >= l)
        index <- l -1
    ## Search for TRUE to the right
    for (i in (index + 1):l) {
        if (sign_change[i] && (i - index) >= min_dist) {
            right_true <- i
            break
        }
    }
    c(left = left_true, right = right_true)
}

#' @title Simple peak detection algorithm on chromatographic data
#'
#' @description
#' 
#' Basic peak picking algorithm that finds a targeted chromatographic peak
#' given based on retention time and intensities and expected retention time
#' of the target signal. The chromatographic signal is expected to contain
#' only signal from a single ion. See details below for information on the
#' algorithm.
#'
#' @details
#'
#' The method:
#'
#' - identifies first all loacal maxima in `intensity`, then
#' - removes local maxima with an intensity lower than 50% of the maximum
#'   intensity in `intensity` and finally
#' - reports the above defined tentative peak apex with an retenton time
#'   closest to the provided `targetRtime`.
#'
#' The closest local minima that are at least 2 data points from the apex
#' position are reported as the left and right margin of the peak.
#' 
#' @param rtime `numeric` with retention times.
#' 
#' @param intensity `numeric` (same length as `rtime`) with the signal
#'     intensities.
#' 
#' @param targetRtime `numeric(1)` with the expected retention time of the
#'     target peak.
#'
#' @param .check `logical(1)` whether input parameters should be checked. Use
#'     `.check = FALSE` only if validity of the input parameters was checked
#'     in an upstream function. 
#' 
#' @author Pablo Vangeenderhuysen
#' 
#' @return
#'
#' @importFrom xcms imputeLinInterpol
#' 
#' @export
find_peak_points <- function(rtime = numeric(), intensity = numeric(),
                             targetRtime = numeric(), .check = TRUE) {
    if (.check) {
        if (!length(intensity))
            stop("'intensity' is of length 0")
        if (length(rtime) != length(intensity))
            stop("'rtime' and 'intensity' are expected to have the same length")
        if (length(targetRtime) != 1L)
            stop("'targetRtime' is expected to be a single numeric")
    }
    ## Compute the derivative of the vector jo: is here something missing?
    if (anyNA(intensity))
        intensity <- imputeLinInterpol(intensity)
    di <- diff(intensity)
    sign_changes <- c(FALSE, diff(di > 0) != 0, FALSE) # same length than ints
    peak_index <- which.max(intensity)
    all_local_max <- which(diff(sign(di)) == -2) + 1
    ## Delete local maxima with an intensity lower than 50% max int
    local_max <- all_local_max[intensity[all_local_max] >
                               0.5 * intensity[peak_index]]
    ## Find the local maximum closest to the targetRtime 
    differences <- abs(rtime[local_max] - targetRtime)
    peak_index <- local_max[which.min(differences)]
    ## Find the left and right border points from the peak
    border <- .find_peak_border(sign_changes, peak_index)

    c(border, peak_index = peak_index)
}


