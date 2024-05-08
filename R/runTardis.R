#' Run T.A.R.D.I.S.
#'
#' Launches the GUI of T.A.R.D.I.S. that allows to input all parameters in an intuitive way to perform targeted peak detection.
#'
#'@import shiny
#'@import shinyFiles
#'@import shinyjs
#'@import jrc
#'@import shinybusy
#'@import bslib
#'@import ggplot2
#'@import stringr
#'@import magrittr
#'@import tibble
#'@import readxl
#'@import openxlsx
#'@import MsExperiment
#'@import BiocParallel
#'@import Spectra
#'@import RColorBrewer
#'@import signal
#'@import xcms
#'@import pracma
#' @export
runTardis <- function() {
  appDir <- system.file("tardis_app", package = "TARDIS")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `TARDIS`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
