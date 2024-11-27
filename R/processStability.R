#' Check Process Stability
#'
#' This function generates an x-bar chart of the variable measured from the reference process and counts the number of out of control points.
#' @section Details:
#' Data must be in the following form to work: it must have n rows (number of observations per subgroup) and m columns (number of subgroups) where the columns are in chronological order of measurement. Each cell must contain a numeric value.
#' @param data values derived from chronological experiment.
#' @return
#'   \describe{
#'     \item{OOC}{The number of out of control points and x-bar chart.}
#'    }
#' @export

processStability = function(data, show.plot = TRUE) {

  xBar = colMeans(data, na.rm=TRUE)
  xBarBar = mean(xBar, na.rm=TRUE)
  R = apply(data, 2, max, na.rm=TRUE)-apply(data, 2, min, na.rm=TRUE)
  RBar = mean(R)
  n = nrow(data)

  A2_Values = c('2' = 1.88, '3' = 1.02, '4' = 0.73, '5' = 0.58, '6' = 0.48, 
                 '7' = 0.42, '8' = 0.37, '9' = 0.34, '10' = 0.31)
  if(as.character(n) %in% names(A2_Values)){
    A2 = A2_Values[as.character(n)]
  } else {
    stop("Error. number of row values is not implemented")
  }

  LCL = xBarBar - A2*RBar
  UCL = xBarBar + A2*RBar

  if(show.plot){
    old_par = par(no.readonly = TRUE)
    par(mar = c(4, 3, 3, 2) + 0.1, xpd = FALSE)
    plt = plot(xBar, xlab = "Time Order", ylab = "X Bar", main = "Control Chart", ylim = c(min(c(xBar, LCL)) - (min(c(xBar, LCL))*.05), max(c(xBar, UCL)) + (max(c(xBar, UCL))*.05)), pch = 19)
    abline(h = LCL, lty = 2, lwd = 2)
    abline(h = UCL, lty = 2, lwd = 2)
    abline(h = xBarBar)
    par(old_par)
  }

  OOC = sum(xBar < LCL | xBar > UCL)

  return(OOC)

}
