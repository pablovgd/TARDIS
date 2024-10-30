#' Find_true_occurence
#'
#' Find peak borders
#'
#' @param signs Vector (TRUE/FALSE) with indices where derivative changes sign
#' (TRUE) or doesn't change sign (FALSE)
#' @param vector Vector with intensities
#' @param point Index of the peak
#'
#' @return
#' @export
#'


find_true_occurrence <- function(signs,vector,point) {
  left_true <- NA
  right_true <- NA

  # If peak is at the leftmost border, then left border is at 1
  #   --> the if in the for loop below will always be FALSE
  if (point == 1) {
    point = point + 1
  }

  # Search for TRUE to the left
  for (i in (point - 1):1) {
    if (signs[i] == TRUE && (point - i) > 2) {
      left_true <- i
      break
    } else {
      left_true <- 1
    }
  }

  # If peak is at the rightmost border, then right border is at vector length
  #   --> the if in the for loop below will always be FALSE
  if (point >= length(signs)) {
    point <- length(signs) - 1
  }

  # Search for TRUE to the right
  for (i in (point + 1):length(signs)) {
    if (signs[i] == TRUE && (i - point) > 2) {
      right_true <- i
      break
    } else {
      right_true <- length(vector)
    }
  }

  result <- list(left_true = left_true, right_true = right_true)
  return(result)
}

#' Find peak points
#'
#' Basic peak picking algorithm that finds a targeted chromatographic peak and
#' its borders given the data and expected retention time of the target.
#'
#' @param rtvector Vector with retention times.
#' @param vector Vector with intensities.
#' @param searchrt Expected retention time of peak.

#' @return
#' @export
#'

find_peak_points <- function(rtvector, vector, searchrt) {
  #Compute the derivative of the vector and find the indices where the
  #  derivative changes sign
  sign_changes <-  c(FALSE, diff(diff(vector)>0)!=0)

  #Find global maximum
  peak_index <- which.max(vector)

  # Find the local maxima
  all_local_max <- which(diff(sign(diff(vector)))==-2)+1

  # Delete local maxima with an intensity lower than 50% of the global maximum
  local_max <- intersect(which(vector > 0.5 * max(vector[all_local_max])),
                         all_local_max)

  # Find the maximum intensity closest to the searchrt, and which is at least
  #   50% of the maximum intensity in the window
  maxrts <- rtvector[c(peak_index, local_max)]
  differences <- abs(maxrts - searchrt)

  # Find the position (index) with the minimum difference
  closest_position <- which.min(differences)
  peak_index <- which(rtvector == maxrts[closest_position])

  # Find the left and right border points from the peak
  occur <- find_true_occurrence(sign_changes, vector, peak_index)

  if (is.na(occur$left_true)) {
    occur$left_true = 1
  }
  if (is.na(occur$right_true)) {
    occur$right_true = length(vector)
  }

  return(list(left = occur$left_true,
              right = occur$right_true,
              foundrt = rtvector[peak_index],
              peakindex = peak_index))
}
