#' Finding Optimal Sampling Strategy
#'
#' This function runs simulations which finds the optimal sampling strategy given various experimental constraints. 
#' @param data dataframe derived from simulateExperiment function 
#' @param wVec a 3 element vector which represents the weights for each desirability function. The first index corresponds to power, the second to sample size, and the third to the false positive rate. Each value must be between 0 and 1 and the sum of all 3 must equal 1. 
#' @return Optimal sampling strategy given by N_1 and T, the final power level, expected sample size (given the two processes are not equivalent), and false positive rate of this design  
#' @export

optimalStrategy = function(data, wVec, L1 = 0, T1 = 1, T2 = 1, U2 = 41, T3 = 0, U3 = 1, k = 2){
  
  dataSampleSize= aggregate(SampleSize~T+n1, data = data, FUN = mean)
  which.min(dataSampleSize$SampleSize)
  
  dataPower = aggregate(Power~T+n1, data = data, FUN = mean)
  dataSampleSize$Power = dataPower$Power
  
  datFN = aggregate(falsePositive~T+n1, data = data, FUN = mean)
  dataSampleSize$falsePositive = datFN$falsePositive
  
  dSampleSize = c()
  dPower = c()
  dFP = c()
  D = c()

  for(i in 1:nrow(dataSampleSize)){
    
    x1 = dataSampleSize$Power[i]
    if(x1<L1){
      dPower[i] = 0
    } else if(x1>T1){
      dPower[i] = 1
    } else {
      dPower[i] = ((x1-L1)/(T1-L1))^k
    }
    
    x2 = dataSampleSize$SampleSize[i]
    if(x2<T2){
      dSampleSize[i] = 1
    } else if(x2>U2){
      dSampleSize[i] = 0
    } else {
      dSampleSize[i] = ((x2-U2)/(T2-U2))^k
    }
    
    x3 = dataSampleSize$falsePositive[i]
    if(x3<T3){
      dFP[i] = 1
    } else if(x3>U3){
      dFP[i] = 0
    } else {
      dFP[i] = ((x3-U3)/(T3-U3))^k
    }
    
    D[i] = dSampleSize[i]^wVec[1] * dPower[i]^wVec[2] * dFP[i]^wVec[3]
    
  }

  index = which.max(D)
  finalPower = dataSampleSize$Power[index]
  finalSampleSize = dataSampleSize$SampleSize[index]
  finalFP = dataSampleSize$falsePositive[index]
  N1 = dataSampleSize$n1[index]
  T = dataSampleSize$T[index]
  
  return(list(finalPower, finalSampleSize, finalFP, N1, T))

}