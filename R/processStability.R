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

  xBar = colMeans(data)
  xBarBar = mean(xBar)
  R = apply(data, 2, max)-apply(data, 2, min)
  RBar = mean(R)
  n = nrow(data)
  if(n == 2){
    A2 = 1.88
  } else if(n == 3){
    A2 = 1.02
  } else if(n == 4){
    A2 = .73
  } else if(n == 5){
    A2 = .58
  } else if(n == 6){
    A2 = .48
  } else if(n == 7){
    A2 = .42
  } else if(n == 8){
    A2 = .37
  } else if(n == 9){
    A2 = .34
  } else if(n == 10){
    A2 = .31
  }

  LCL = xBarBar - A2*RBar
  UCL = xBarBar + A2*RBar

  if(show.plot){
    plt = plot(xBar, xlab = "Time Order", ylab = "X Bar", main = "Control Chart", ylim = c(min(c(xBar, LCL)) - (min(c(xBar, LCL))*.05), max(c(xBar, UCL)) + (max(c(xBar, UCL))*.05)), pch = 19)
    abline(h = LCL, lty = 2, lwd = 2)
    abline(h = UCL, lty = 2, lwd = 2)
    abline(h = xBarBar)
  }

  OOC = sum(xBar < LCL | xBar > UCL)

  return(OOC)

}
