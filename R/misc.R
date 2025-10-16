#' Helper function to save parameter inputs as .csv by collapsing strings
#'
#' @param x `
#'
#'
#' @return collapsed strings
#'
#' @author Pablo Vangeenderhuysen
#'
#' @noRd
.collapse_safe <- function(x, collapse = ",", null_string = "") {
  if (is.null(x)) {
    null_string
  } else if (length(x) == 0) {
    ""
  } else {
    paste(x, collapse = collapse)
  }
}
