#' Finding Optimal Sampling Strategy
#'
#' Finds the optimal sampling strategy given results from simulations.
#' @section Details:
#' There are 3 metrics that can be considered when selecting optimal sampling strategies; power, the false positive rate, and the expected sample size (given that the reference and candidate processes are not equivalent). The strategy which best trades-off these metrics can be considered the optimal sampling strategy. This function performs multi-objective optimization with desirability functions, converting each metric into a desirability value and then finding the sampling strategy that maximizes the product of the desirabilities taken to the weights associated with each desirability. The metrics are derived from simulations given by the simulateExperiment function.
#' @param data dataframe containing numerical results from simulateExperiment function.
#' @param wVec a 3 element numeric vector which represents the weights for each desirability function. The first index corresponds to power, the second to sample size, and the third to the false positive rate. Each value must be between 0 and 1 and the sum of all 3 must equal 1.
#' @param L1 a numeric value representing the lowest acceptable level of power.
#' @param T1 a numeric value representing the target level of power.
#' @param T2 a numeric value representing the target sample size.
#' @param U2 a numeric value representing the maximum acceptable sample size.
#' @param T3 a numeric value representing the target false positive rate.
#' @param U3 a numeric value representing the maximal acceptable false positive rate.
#' @param k a numeric value giving the concavity of each desirability function.
#'
#' @return A list containing the metrics and experimental parameters of the optimal sampling strategy.
#'   \describe{
#'     \item{finalPower}{The value of the optimal power level.}
#'     \item{finalSampleSize}{The value of the optimal expected sample size.}
#'     \item{finalFP}{The value of the optimal false positive rate.}
#'     \item{N1}{The value of the optimal initial sample size.}
#'     \item{T}{The value of the optimal number of sequential samples.}
#'   }
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
