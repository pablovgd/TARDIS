ms_data_frame <- function(x) {
  pks <- peaksData(x)
  npks <- vapply(pks, nrow, integer(1))
  res <- as.data.frame(do.call(rbind, pks))
  res$rtime <- rep(rtime(x), npks)
  res
}