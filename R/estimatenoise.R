
estimateBaselineNoise = function(eic, peakwidth) {
  #Accepts an EIC with columns "i" and "rt" and adds columns "i.sg", "baseline.region", "baseline", "noise.sd"
  # This function wraps up baseline estimation and noise estimation.  Baseline estimation is performed on SG smoothed EIC traces.  Regions with a slope greater than minslope.peak are excluded from the baseline estimation.  The smoothed, nonpeak EIC is interpolated and taken as the baseline.  The original, nonpeak EIC is applied a rolling sd calculation and that is taken as the SD of the noise at each point.
  
  # min.peakwidth is used to set the windows for smoothing and noise SD estimation
  
  eicr = matrix(numeric(), nrow = nrow(eic), ncol = ncol(eic) + 3)
  colnames(eicr) = c(colnames(eic), "i.sg", "baseline", "noise.sd")
  eicr[,1:ncol(eic)] = unlist(eic)
  
  # Turn supplied peakwidth into scans
  scantimes = diff(eic[,"rt"])
  if (diff(quantile(scantimes, c(0.05, 0.95)))/mean(scantimes) > 1) {
    warning("Scan rate varies by more than 1%.  Could cause errors in peak detection. Peak width used will be incorrect.")
  }
  scanwidth = round(peakwidth / mean(scantimes))
  
  # Compute a smoothed EIC
  sg.window = round(scanwidth[1]/2) %>% { if (.%%2 == 0) .+1 else . }
  eicr[,"i.sg"] = signal::sgolayfilt(eicr[,"i"], p = 2, n = sg.window)
  
  # Find the SD of the residuals of the smoothing => Noise SD
  tmpeic = eicr[,"i"] %>% {.[.==0] = eicr[.==0,"i.sg"]*0.5; .}
  rollav = stats::filter(eicr[,"i"], 1/rep(sg.window, sg.window)) %>% { .[is.na(.)] = eicr[,"i"][is.na(.)]; . }
  
  #eic$noise.local.sd = zoo::rollapply(tmpeic - ( rollav %>% {.[. <= 0] = 0; .} ), sg.window, "sd", fill=NA) %>% { .[is.na(.)] = 0; . }
  foo = tmpeic - rollav
  eicr[,"noise.sd"] = TTR::runSD(foo, sg.window) %>% { .[is.na(.)] = 0; . }
  #eic$noise.local.sd = zoo::rollmean(eic$noise.local.sd, round(sg.window/2), fill = NA)
  
  eicr[,"baseline"] = baseline::baseline(matrix(eicr[,"i"], nrow = 1), method="als", p = 0.1, lambda = ceiling(scanwidth[1]/4))@baseline
  
  eicr  # returns a matrix with columns "i", "i.sg", "rt", "baseline.region", "baseline", "noise.baseline.sd", "noise.local.sd"
}