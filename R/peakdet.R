find_true_occurrence <- function(signs,vector,point) {
  left_true <- NA
  right_true <- NA
  
  if(point == 1){
    point = point + 1
  }
  
 

  # Search for TRUE to the left --> if sign changes and at least doens't change at the next point 
  for (i in (point - 1):1) {
    if (signs[i] == TRUE &&  point - i > 2) {
      left_true <- i
      break
    } else{ left_true <- 1}
  }
  
  if(point >= length(signs)){
    point = length(signs) -1
  }
  
  # Search for TRUE to the right
  for (i in (point + 1):length(signs)) {
    if (signs[i] == TRUE && i - point > 2) {
      right_true <- i
      break
    } else{ right_true <- length(vector)}
  }
  
  result <- list(left_true = left_true, right_true = right_true)
  return(result)
}

find_peak_points <- function(rtvector, vector, searchrt) {
  # Compute the derivative of the vector

  
  # Find the indices where the derivative changes sign
  sign_changes <-  c(FALSE, diff(diff(vector)>0)!=0)
  
  # Find absolute maximum
  peak_index <- which.max(vector)
  
  # Find the local maxima
  
  all_local_max <- which(diff(sign(diff(vector)))==-2)+1
  
  
  # Delete local maxima with an intensity lower than 50k
  
  local_max <- c()
  
  for(lmax in all_local_max){
    if(vector[lmax] > 0.5* max(vector[all_local_max])){
      local_max <- cbind(local_max,lmax)
    }
  }
  
  
  # Find the max intensity closest to the searchrt which is at least 50% of the max intensity in the window
  
  maxrts <- rtvector[c(peak_index,local_max)]
  
 
  differences <- abs(maxrts - searchrt)
  
  # Find the position (index) with the minimum difference
  closest_position <- which.min(differences)
  
  peak_index <- which(rtvector == maxrts[closest_position])
  
  # Find the left and right border points from the peak
  occur <- find_true_occurrence(sign_changes,vector,peak_index)
  
  if(is.na(occur$left_true)){
    occur$left_true = 1
  }
  if(is.na(occur$right_true)){
    occur$right_true = length(vector)
  }
  
  return(list(left = occur$left_true, right = occur$right_true,foundrt = rtvector[peak_index],peakindex = peak_index))
}
