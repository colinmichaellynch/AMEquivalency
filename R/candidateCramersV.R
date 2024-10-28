#' Candidate Cramer's V
#'
#' Calculates the normalized Cramer's V of the candidate process based on the limits of the reference bins.
#' @param chiSquare a numeric value giving the goodness of fit of the candidate samples on the reference bins
#' @param b a numeric value indicating the number of bins
#' @return
#'   \describe{
#'     \item{N}{A numeric value indicating the sample size.}
#'    }
#' @export

candidateCramersV = function(chiSquare, b, N){

  return(sqrt(chiSquare/N)/sqrt(b-1))

}
