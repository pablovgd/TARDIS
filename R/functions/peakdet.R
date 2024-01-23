find_true_occurrence <- function(vector, point) {
  left_true <- NA
  right_true <- NA
  
  if(point == 1){
    point = point + 1
  }
  
 

  # Search for TRUE to the left
  for (i in (point - 1):1) {
    if (vector[i] == TRUE) {
      left_true <- i
      break
    }
  }
  
  if(point >= length(vector)){
    point = length(vector) -1
  }
  
  # Search for TRUE to the right
  for (i in (point + 1):length(vector)) {
    if (vector[i] == TRUE) {
      right_true <- i
      break
    }
  }
  
  result <- list(left_true = left_true, right_true = right_true)
  return(result)
}

find_peak_points <- function(vector) {
  # Compute the derivative of the vector

  
  # Find the indices where the derivative changes sign
  sign_changes <-  c(FALSE, diff(diff(vector)>0)!=0)
  
  # Find the index of the peak
  peak_index <- which.max(vector)
  
  # Find the left and right points from the peak
  occur <- find_true_occurrence(sign_changes,peak_index)
  
  if(is.na(occur$left_true)){
    occur$left_true = 1
  }
  if(is.na(occur$right_true)){
    occur$right_true = length(vector)
  }
  
  return(list(left = occur$left_true, right = occur$right_true))
}
