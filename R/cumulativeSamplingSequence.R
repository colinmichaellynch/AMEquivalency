#' Cumulative Sampling Sequence
#'
#' This function calculates the sample size of sequential samples for the equivalency experiment. Sample sizes are given in terms of cumulative sample size. 
#' @param N_1 initial sample size
#' @param T number of sequential samples
#' @return N_T sampling sequence
#' @export

cumulativeSamplingSequence = function(N_1, T) {
  
  t = 1:T
  N_T = round((2*t*T*N_1+t^2-t^2*N_1 + t - t*N_1 - 2 + 2*N_1)/(2*T))
  return(N_T)
  
}