#' Candidate Process Chi Squared Value 
#'
#' This function calculates the chi squared value of the candidate process as it fits to the bins defined by the reference process
#' Data must be in the following form to work: it must have n rows (number of observations per subgroup) and m columns (number of subgroups) where the columns are in chronological order of measurement. Each cell must contain a numeric value. bins must be defined by quantile function (see tutorial)
#' @param data values derived from chronological experiment
#' @param bins percentiles of reference distribution
#' @return chi squared value
#' @export

candidateChiSquared = function(data, bins, show.plot = TRUE){
  
  vec = c(as.matrix(data))
  N = length(vec)
  dataBinned = cut(vec, breaks = bins, include.lowest = TRUE)
  observed = table(dataBinned)
  chiSquared = sum((observed-(N/b))^2/(N/b))
  
  if(show.plot){
    dataPlot = data.frame(Bins = names(observed), Count = as.numeric(observed))
    barplot(height = dataPlot$Count, names.arg = dataPlot$Bins, xlab = "", ylab = "Count", las = 2)  
  }
  
  return(chiSquared)
  
}
  