#' Candidate Process Chi Squared Value
#'
#' This function calculates the chi squared value of the candidate process as it fits to the bins defined by the reference process.

#' @section Details:
#' Data must be in the following form: it must have n rows (number of observations per subgroup) and m columns (number of subgroups) where the columns are in chronological order of measurement. Each cell must contain a numeric value. Bins must be defined by quantile function (see tutorial).
#' @param data a numeric vector or matrix containing values derived from candidate process measurements.
#' @param bins a numeric vector containing percentiles of reference distribution.
#' @param show.plot a logical indicating whether to show the x-bar chart.
#' @return
#'   \describe{
#'     \item{chiSquared}{A value representing the candidate chi square value}
#'    }
#' @export

candidateChiSquared = function(data, bins, show.plot = TRUE){

  vec = c(as.matrix(data))
  N = length(vec)
  dataBinned = cut(vec, breaks = bins, include.lowest = TRUE)
  observed = table(dataBinned)
  chiSquared = sum((observed-(N/b))^2/(N/b))

  if(show.plot){
    old_par = par(no.readonly = TRUE)
    par(mar = c(6, 5, 3, 1) - 0.25)
    dataPlot = data.frame(Bins = names(observed), Count = as.numeric(observed))
    barplot(height = dataPlot$Count,names.arg = dataPlot$Bins,xlab = "",ylab = "Count",las = 2,cex.names = 0.8)
    par(old_par)
  }

  return(chiSquared)

}
