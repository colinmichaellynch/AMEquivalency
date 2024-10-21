#' Truncated Geometric Distribution 
#'
#' This function generates the PMF of a truncated geometric distribution
#' @param bins a vector containing the indeces for each bin 
#' @param p rho governing how much of the data is either uniformly distributed or clustered to one side of the distribution
#' @param b the number of bins
#' @return vector containing PMF values for 1 to b bins corresponding to truncated geometric distribution set by p
#' @export

truncatedGeometricPMF = function(bins, p, b) {
  
  if (bins > b) return(0)
  pmfUntruncated = (1 - p)^(bins - 1) * p
  normalizationFactor = sum((1 - p)^(0:(b - 1)) * p)
  pmfTruncated = pmfUntruncated / normalizationFactor
  return(pmfTruncated)
  
}