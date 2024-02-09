chrom_peak_intensity_centWave <- function(x, rt, peakArea,
                                           mzCenterFun = "weighted.mean",
                                           sampleIndex = integer(),
                                           cn = character(), ...) {
  cn = c(cn,"sn")
  res <- matrix(NA_real_, ncol = length(cn), nrow = nrow(peakArea))
  rownames(res) <- rownames(peakArea)
  colnames(res) <- cn
  res[, "sample"] <- sampleIndex
  res[, c("rtmin", "rtmax", "mzmin", "mzmax")] <-
    peakArea[, c("rtmin", "rtmax", "mzmin", "mzmax")]
  for (i in seq_len(nrow(res))) {
    rtr <- peakArea[i, c("rtmin", "rtmax")]
    keep <- which(MsCoreUtils::between(rt, rtr))
    if (length(keep)) {
      xsub <- lapply(x[keep], xcms:::.pmat_filter_mz,
                     mzr = peakArea[i, c("mzmin", "mzmax")])
      ## length of xsub is the number of spectra, the number of peaks can
      ## however be 0 if no peak was found. Maybe we should/need to
      ## consider adding 0 or NA intensity for those.
      mat <- do.call(rbind, xsub)
      if (nrow(mat)) {
        ## can have 0, 1 or x values per rt; repeat rt accordingly
        rts <- rep(rt[keep], vapply(xsub, nrow, integer(1L)))
        maxi <- which.max(mat[, 2L])[1L]
        mmz <- do.call(mzCenterFun, list(mat[, 1L], mat[, 2L]))
        if (is.na(mmz)) mmz <- mat[maxi, 1L]
        res[i, c("rt", "mz", "maxo", "into","sn")] <- c(
          rts[maxi], mmz, mat[maxi, 2L],
          sum(mat[, 2L], na.rm = TRUE) *
            ((rtr[2L] - rtr[1L]) / max(1L, (length(keep) - 1L))),
          mat[maxi, 2L]/ xcms:::estimateChromNoise(mat[, 2L])
        )
        if ("beta_cor" %in% cn)
          res[i, c("beta_cor", "beta_snr")] <- xcms:::.get_beta_values(
            mat[, 2L], rts)
      }
    }
  }
  res[!is.na(res[, "maxo"]), , drop = FALSE]
}

xmse_integrate_chrom_peaks <- function(x, pal, msLevel = 1L,
                                        intFun = chrom_peak_intensity_centWave,
                                        mzCenterFun = xcms:::mzCenter.wMean,
                                        param = MatchedFilterParam(),
                                        BPPARAM = bpparam(), ...) {
  keep <- which(msLevel(spectra(x)) == msLevel)
  f <- as.factor(fromFile(x)[keep])
  if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
  else rt <- rtime(spectra(x))[keep]
  cn <- colnames(chromPeaks(x))
  res <- bpmapply(split(Spectra::peaksData(filterMsLevel(spectra(x), msLevel),
                                  f = factor()), f),
                  split(rt, f),
                  pal,
                  as.integer(names(pal)),
                  FUN = intFun,
                  MoreArgs = list(mzCenterFun = mzCenterFun, cn = cn,
                                  param = param),
                  SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
  do.call(rbind, res)
}



manualChromPeaks_2 = function (object, chromPeaks = matrix(numeric()), 
                                           samples = seq_along(object), msLevel = 1L, chunkSize = 2L, 
                                           BPPARAM = bpparam()) 
{
  if (length(msLevel) > 1L) 
    stop("Can only add peaks from one MS level at a time.")
  if (is.data.frame(chromPeaks)) 
    chromPeaks <- as.matrix(chromPeaks)
  if (!nrow(chromPeaks)) 
    return(object)
  if (!all(c("mzmin", "mzmax", "rtmin", "rtmax") %in% 
           colnames(chromPeaks))) 
    stop("'chromPeaks' lacks one or more of the required colums ", 
         "\"mzmin\", \"mzmax\", \"rtmin\" and \"rtmax\".")
  chromPeaks <- chromPeaks[, c("mzmin", "mzmax", "rtmin", 
                               "rtmax")]
  if (!all(samples %in% seq_along(object))) 
    stop("'samples' out of bounds")
  if (hasFeatures(object)) 
    object <- dropFeatureDefinitions(object)
  pal <- lapply(samples, function(z) chromPeaks)
  names(pal) <- samples
  chunks <- split(samples, ceiling(seq_along(samples)/chunkSize))
  pb <- progress_bar$new(format = paste0("[:bar] :current/:", 
                                         "total (:percent) in ", ":elapsed"), total = length(chunks) + 
                           1L, clear = FALSE)
  pb$tick(0)
  res <- lapply(chunks, function(z, ...) {
    pb$tick()
    xmse_integrate_chrom_peaks(xcms:::.subset_xcms_experiment(object, 
                                                        i = z, keepAdjustedRtime = TRUE, ignoreHistory = TRUE), 
                                pal = pal[as.character(z)], msLevel = msLevel, 
                                BPPARAM = BPPARAM)
  })
  res <- do.call(rbind, res)
  nr <- nrow(res)
  maxi <- max(0, as.integer(sub("CP", "", rownames(xcms:::.chromPeaks(object)))))
  rownames(res) <- xcms:::.featureIDs(nr, "CP", maxi + 1)
  pkd <- data.frame(ms_level = rep(as.integer(msLevel), 
                                   nr), is_filled = rep(FALSE, nr))
  rownames(pkd) <- rownames(res)
  pb$tick()
  object@chromPeakData <- MsCoreUtils:::rbindFill(object@chromPeakData, 
                                    pkd)
  object@chromPeaks <- MsCoreUtils:::rbindFill(object@chromPeaks, res)
  validObject(object)
  object
}




