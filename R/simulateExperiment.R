#' Experimental Simulations
#'
#' This function simulates experiments to determine the power, expected sample size, and false positive rates of various sampling strategies.
#' @section Details:
#' Simulates experiments where the null hypothesis is either true or not true. The null hypothesis is that the candidate and reference processes are equivalent. The probability that an experiment detects a true difference between these two processes is the power of that design. The probability that the experiment indicates that the two processes are not equivalent when in fact they are is the false positive rate. The expected sample size is the average sample size of the experiment given that the two processes are not equivalent (the experiment stops early because there is enough evidence to suggest that the null hypothesis is not true.) Experimental designs are set by the function possibleSamplingStrats.
#' @param b the number of bins used to characterize the reference distribution
#' @param alpha the significance level
#' @param n_max the highest possible sample size
#' @param simulations the number values of rho tested in simulations to find optimal sampling strategy
#' @param iterations the number of simulations per parameter value rho
#' @param TVecAcceptable vector containing values of T for possible experiments
#' @param n1VecAcceptable vector containing values of N_1 for possible experiments
#' @return
#'   \describe{
#'     \item{data}{A dataframe containing the power levels, sample sizes, and false positive rates of various combinations of N_1 and T.}
#'    }
#' @export

simulateExperiment = function(b = 8, alpha = .01, n_max = 499, simulations = 100, iterations = 100, TVecAcceptable, n1VecAcceptable){

  bins = 1:b
  chiSquareU = qchisq(p = alpha, df = b-1, lower.tail = FALSE)
  n = 1:n_max
  vTilde = sqrt(n)/n
  A = (sqrt(2)*gamma(b/2))/(sqrt(b-1)*gamma((b-1)/2))

  pVec = c()
  n1Vec = c()
  TVec = c()
  tPass = c()
  powerVec = c()
  falsePositive = c()
  cramersVVec = c()
  sampleSize = c()
  counter = 0

  for(m in 1:length(TVecAcceptable)){

    T = TVecAcceptable[m]
    n_1 = n1VecAcceptable[m]
    s = (1-n_1)/T
    n_t = round(((1:T)*n_1)+((((1:T)-1)*((1:T)+2))/2)*s)

    for(i in 1:simulations){

      counter = counter+1
      p = runif(1, 0, 1)
      probabilities = sapply(bins, truncatedGeometricPMF, p = p, b = b)
      truePosCounter = 0
      TPassTotal = c()
      cramersVPass = c()
      cramersV = c()
      sampleSizeVec = c()

      for(j in 1:iterations){

        for(t in 1:T){

          if(t == 1){
            sampleRaw = sample(1:b, n_t[1], replace = TRUE, prob = probabilities)
          } else {
            sampleRaw = c(sampleRaw, sample(1:b, n_t[t]-n_t[t-1], replace = TRUE, prob = probabilities))
          }

          sample = table(factor(sampleRaw, levels = 1:b))
          chiSquare = sum(((sample-rep(n_t[t]/b, b))^2)/rep(n_t[t]/b, b))
          cramersV[t] = sqrt(chiSquare/n_t[t])/sqrt(b-1)

        }

        cramersVUpper = sqrt(chiSquareU/n_t)/sqrt(b-1)

        if(any(cramersV>cramersVUpper)){
          truePosCounter = truePosCounter + 1
          TPassTotal[j] = min(which(cramersV>cramersVUpper))
          cramersVPass[j] = cramersV[TPassTotal[j]]
          sampleSizeVec[j] = n_t[TPassTotal[j]]
        } else {
          TPassTotal[j] = NaN
          cramersVPass[j] = NaN
          sampleSizeVec[j] = NaN
        }

      }

      pVec[counter] = p
      n1Vec[counter] = n_1
      TVec[counter] = T
      tPass[counter] = median(TPassTotal, na.rm = TRUE)
      cramersVVec[counter] = median(cramersVPass, na.rm = TRUE)
      sampleSize[counter] = median(sampleSizeVec, na.rm = TRUE)
      powerVec[counter] = truePosCounter/iterations
      falsePositive[counter] =NaN

    }

    for(i in 1:simulations){

      counter = counter+1
      probabilities = rep(1/b, b)
      falsePositiveCounter = 0
      TPassTotal = c()
      cramersVPass = c()
      cramersV = c()
      sampleSizeVec = c()

      for(j in 1:iterations){

        for(t in 1:T){

          if(t == 1){
            sampleRaw = sample(1:b, n_t[1], replace = TRUE, prob = probabilities)
          } else {
            sampleRaw = c(sampleRaw, sample(1:b, n_t[t]-n_t[t-1], replace = TRUE, prob = probabilities))
          }

          sample = table(factor(sampleRaw, levels = 1:b))
          chiSquare = sum(((sample-rep(n_t[t]/b, b))^2)/rep(n_t[t]/b, b))
          cramersV[t] = sqrt(chiSquare/n_t[t])/sqrt(b-1)

        }

        cramersVUpper = sqrt(chiSquareU/n_t)/sqrt(b-1)

        if(any(cramersV>cramersVUpper)){
          falsePositiveCounter = falsePositiveCounter + 1
          TPassTotal[j] = min(which(cramersV>cramersVUpper))
          cramersVPass[j] = cramersV[TPassTotal[j]]
          sampleSizeVec[j] = n_t[TPassTotal[j]]
        } else {
          TPassTotal[j] = NaN
          cramersVPass[j] = NaN
          sampleSizeVec[j] = NaN
        }

      }

      pVec[counter] = p
      n1Vec[counter] = n_1
      TVec[counter] = T
      tPass[counter] = NaN
      cramersVVec[counter] = median(cramersVPass, na.rm = TRUE)
      sampleSize[counter] = NaN
      powerVec[counter] = NaN
      falsePositive[counter] = falsePositiveCounter/iterations

    }

    print(paste("Percent complete: ", m/length(TVecAcceptable)*100))

  }

  data = (data.frame(p = pVec, n1 = n1Vec, T = TVec, tPass = tPass, cramersV = cramersVVec, Power = powerVec, SampleSize = sampleSize, falsePositive = falsePositive))

  return(data)

}
