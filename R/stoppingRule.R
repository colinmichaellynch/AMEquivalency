#' When to Stop Sampling
#'
#' This function takes the cramers v of a sample and other aspects of the design to determine whether or not the statistic is outside the limits expected for that design given the null hypothesis that the candidate and reference processes are equivalent.
#' @param cramersV measurment of how closely the candidate process resembles the reference process
#' @param N sample size corresponding to this value of cramers's v
#' @param b the number of bins
#' @param alpha significance level
#' @return
#'   \describe{
#'     \item{stop}{A logical representing the decision to continue sampling. If stop == TRUE, keep going unless t = T. If stop == FALSE, stop sampling.}
#'    }
#' @export

stoppingRule = function(cramersV, N, b, alpha = .01){

  criticalChiSquare = qchisq(alpha, df = N-1)
  upperVal = 1/sqrt(N*(b-1)) * sqrt(criticalChiSquare)
  if(cramersV>upperVal){
    stop = TRUE
  } else {
    stop = FALSE
  }

  return(stop)

}
