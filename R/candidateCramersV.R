#' @export

candidateCramersV = function(chiSquare, b, N){
  
  return(sqrt(chiSquare/N)/sqrt(b-1))
  
}