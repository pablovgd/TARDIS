calculateZigZagIndex <- function(pts){
  
  intPts <- pts
  
  if(length(intPts)>3){
    eic <- intPts
    
    end <- length(eic)
    EPI=max(eic)-(eic[1]+eic[2]+eic[end]+eic[end-1])/4.0;
    
    zig_zag_sum=0.0
    for(i in 2:(end-1)){
      local_zig_zag=(2*eic[i]-eic[i-1]-eic[i+1])^2.0
      zig_zag_sum=zig_zag_sum+local_zig_zag
    }
    
    zig_zag_index = zig_zag_sum/(EPI^2.0*end)
  }else{
    zig_zag_index <- NA
  }
  return(zig_zag_index)
}

calculateGaussianSimilarity <- function(rt, int, rtdb){
 

  if(length(rt) > 2){
    num_peak_pts <- length(rt)
    td <- rt
    d <- int
    mu <- rtdb
    sigma <- tail(rt,1) - rt[1]
    h <- max(int)
    
    fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)
    
    if(class(fit) != "try-error"){
      gaussPts <- as.matrix(fitted(fit))
      gaussPts_std <- (gaussPts-mean(gaussPts))/sd(gaussPts)
      gaussPts_scale <- gaussPts_std/norm(gaussPts_std, type="F")
      
      d <- as.matrix(d)
      peak_intensity_std <- (d-mean(d))/sd(d)
      peak_intensity_scale <- peak_intensity_std/norm(peak_intensity_std, type="F")
      
      gauss_similarity <- sum(gaussPts_scale*peak_intensity_scale)
      
    }else{
      gauss_similarity <- NA
    }
  }else{
    gauss_similarity <- NA
  }
  return(gauss_similarity)
}


qscoreCalculator <- function(rt, int){
  #Check for bogus EICs
  if(length(rt)<5){
    return(list(SNR=NA, peak_cor=NA))
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (rt-min(rt))/(max(rt)-min(rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, int)
  
  
  #Calculate the normalized residuals
  residuals <- int/max(int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(int)-min(int))/sd(norm_residuals*max(int))
  #Return the quality score
  return(list(SNR=SNR, peak_cor=peak_cor))
}
