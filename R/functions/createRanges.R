createRanges <- function(msData,dbData,ppm,rtdev){



### Code to create mz & rt range for compounds ###
  
  mzmed <- dbData$`m/z`
  rtmed <- dbData$tr
  
  deltaTR <- rtdev
  
  mzdeltas <- sapply(mzmed, function(mzmed)
    mzmed * ppm / 10 ^ 6)
  
  #calculate mzrange based on delta mz
  mzRanges <-
    cbind(as.numeric(mzmed) - mzdeltas,
          as.numeric(mzmed) + mzdeltas)
  
  spectra = msData@spectra
  
  #if an upper m/z boundary is lower than the minimum m/z range, set it to the minimum m/z
  indexTemp <-
    which(mzRanges[, 2] < min(spectra@backend@spectraData@listData$basePeakMZ))
  mzRanges[indexTemp, ] <-
    min(spectra@backend@spectraData@listData$basePeakMZ)
  
  #if a lower m/z boundary is higher than the max m/z, set it to the max m/z
  indexTemp <-
    which(mzRanges[, 1] > max(spectra@backend@spectraData@listData$basePeakMZ))
  mzRanges[indexTemp, ] <-
    max(spectra@backend@spectraData@listData$basePeakMZ)
  
  #if an upper limit is higher than the max mz & the lower limit is smaller than the max m/z, set the upper limit to the max m/z
  mzRanges[which(
    mzRanges[, 2] > max(spectra@backend@spectraData@listData$basePeakMZ) &
      mzRanges[, 1] < max(spectra@backend@spectraData@listData$basePeakMZ)
  ), 2] <-
    max(spectra@backend@spectraData@listData$basePeakMZ)
  
  #if a upper limit is larger than the minimum & the lower limit is lower than the minimum, set the lower limit to the min m/z
  mzRanges[which(
    mzRanges[, 2] > min(spectra@backend@spectraData@listData$basePeakMZ) &
      mzRanges[, 1] < min(spectra@backend@spectraData@listData$basePeakMZ)
  ), 1] <-
    min(spectra@backend@spectraData@listData$basePeakMZ)
  
  
  
  rtRanges <-
    cbind(as.numeric(rtmed) - deltaTR / 2,
          as.numeric(rtmed) + deltaTR / 2)
  
  #for al upper RT limits -2 lower than minimum RT, change to lower limit of min rt & upper limit to min + 10
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
  
  
  #for al lower RT limits +2 higher than maximum RT, change lower limit to max -10 and upper limit to max
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
  
  
  #for upper limits higher than max and lower limit lower than max, change upper limit to max
  rtRanges[which(
    rtRanges[, 2] > max(spectra@backend@spectraData@listData$rtime) &
      rtRanges[, 1] < max(spectra@backend@spectraData@listData$rtime)
  ), 2] <- max(spectra@backend@spectraData@listData$rtime)
  #for upper limits higher than min and lower limit lower than min, change lower limit to min
  rtRanges[which(
    rtRanges[, 2] > min(spectra@backend@spectraData@listData$rtime) &
      rtRanges[, 1] < min(spectra@backend@spectraData@listData$rtime)
  ), 1] <- min(spectra@backend@spectraData@listData$rtime)

  return(list(mzRanges,rtRanges))
  
}