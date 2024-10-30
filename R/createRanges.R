#' Create m/z and retention time ranges for target compounds
#'
#' Creates ranges around given m/z and retention time based on given data and
#' allowed deviance.
#'
#' @param msData [MsExperiment()] object
#' @param dbData Target database, output of [createTargetList()]
#' @param ppm Allowed deviance in ppm around given m/z value
#' @param rtdev Allowed deviance in seconds of retention time. Defines the
#' search window in the time dimension.
#'
#' @return A list containing the m/z and retention time ranges for all given
#' target compounds
#' @export
#'


createRanges <- function(msData, dbData, ppm, rtdev){

  mzmed <- dbData$`m/z`
  rtmed <- dbData$tr

  deltaTR <- rtdev

  mzdeltas <- sapply(mzmed, function(mzmed)
    mzmed * ppm / 10 ^ 6)

  #Calculate mzrange based on delta mz
  mzRanges <-
    cbind(as.numeric(mzmed) - mzdeltas,
          as.numeric(mzmed) + mzdeltas)

  spectra <- msData@spectra

  #If an upper m/z boundary is lower than the minimum m/z range, set it to the
  #  minimum m/z
  indexTemp <-
    which(mzRanges[, 2] < min(spectra@backend@spectraData@listData$basePeakMZ))
  mzRanges[indexTemp, ] <-
    min(spectra@backend@spectraData@listData$basePeakMZ)

  #If a lower m/z boundary is higher than the max m/z, set it to the max m/z
  indexTemp <-
    which(mzRanges[, 1] > max(spectra@backend@spectraData@listData$basePeakMZ))
  mzRanges[indexTemp, ] <-
    max(spectra@backend@spectraData@listData$basePeakMZ)

  #If an upper limit is higher than the max m/z & the lower limit is smaller
  #  than the max m/z, set the upper limit to the max m/z
  mzRanges[which(
    mzRanges[, 2] > max(spectra@backend@spectraData@listData$basePeakMZ) &
      mzRanges[, 1] < max(spectra@backend@spectraData@listData$basePeakMZ)
  ), 2] <-
    max(spectra@backend@spectraData@listData$basePeakMZ)

  #If an upper limit is larger than the minimum m/z & the lower limit is lower
  #  than the minimum, set the lower limit to the min m/z
  mzRanges[which(
    mzRanges[, 2] > min(spectra@backend@spectraData@listData$basePeakMZ) &
      mzRanges[, 1] < min(spectra@backend@spectraData@listData$basePeakMZ)
  ), 1] <-
    min(spectra@backend@spectraData@listData$basePeakMZ)



  rtRanges <-
    cbind(as.numeric(rtmed) - deltaTR / 2,
          as.numeric(rtmed) + deltaTR / 2)

  #For all upper RT limits -2 lower than minimum RT, change to lower limit of
  #  min rt & upper limit to min + 10
  indexTemp <-
    which(rtRanges[, 2] - 2 < min(spectra@backend@spectraData@listData$rtime))
  rtRanges[indexTemp, ] <-
    cbind(rep(
      min(spectra@backend@spectraData@listData$rtime),
      length(indexTemp)
    ), rep(
      min(spectra@backend@spectraData@listData$rtime) + 10,
      length(indexTemp)
    ))


  #For all lower RT limits +2 higher than maximum RT, change lower limit to max
  #  -10 and upper limit to max
  indexTemp <-
    which(rtRanges[, 1] + 2 > max(spectra@backend@spectraData@listData$rtime))
  rtRanges[indexTemp, ] <-
    cbind(rep(
      max(spectra@backend@spectraData@listData$rtime) - 10,
      length(indexTemp)
    ), rep(
      max(spectra@backend@spectraData@listData$rtime),
      length(indexTemp)
    ))


  #For upper limits higher than max and lower limit lower than max, change upper
  #  limit to max
  rtRanges[which(
    rtRanges[, 2] > max(spectra@backend@spectraData@listData$rtime) &
      rtRanges[, 1] < max(spectra@backend@spectraData@listData$rtime)
  ), 2] <- max(spectra@backend@spectraData@listData$rtime)
  #For upper limits higher than min and lower limit lower than min, change lower
  #  limit to min
  rtRanges[which(
    rtRanges[, 2] > min(spectra@backend@spectraData@listData$rtime) &
      rtRanges[, 1] < min(spectra@backend@spectraData@listData$rtime)
  ), 1] <- min(spectra@backend@spectraData@listData$rtime)

  return(list(mzRanges,rtRanges))

}
