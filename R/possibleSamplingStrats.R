#' Set of Possible Designs
#'
#' This function finds all of the possible sampling strategies whose final sample size is equal to a constant set by the user (N_1T_Max). 
#' @param N_1T_Max final sample size of experiment given the reference and candidate processes are equivalent
#' @return List of N_1 and T values which produce strategies whose final sample size is equal to N_1T_Max
#' @export

possibleSamplingStrats = function(N_1T_Max=40){
 
  TVecAcceptable = c()
  n1VecAcceptable = c()
  counter = 0
  for(n_1 in 1:N_1T_Max){
    
    for(T in 1:40){
      
      s = (1-n_1)/T
      n_t = round(((1:T)*n_1)+((((1:T)-1)*((1:T)+2))/2)*s)
      
      if(n_t[T] == N_1T_Max){
        
        counter = counter + 1
        TVecAcceptable[counter] = T
        n1VecAcceptable[counter] = n_1
        
      }
      
    }
    
  }
  
  return(list(TVecAcceptable, n1VecAcceptable))
   
}