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
simulateExperiment <- function(b = 8, alpha = 0.01, n_max = 499, simulations = 100, iterations = 100, TVecAcceptable, n1VecAcceptable) {
  
  # Set constants
  chiSquareU <- qchisq(p = alpha, df = b - 1, lower.tail = FALSE)
  bins <- 1:b
  A <- (sqrt(2) * gamma(b / 2)) / (sqrt(b - 1) * gamma((b - 1) / 2))
  
  # Helper function: Simulate one experiment
  # Cleans variable workspace after experiment simulation complete
  simulate_one_experiment <- function(probabilities, T, n_t, chiSquareU, iterations, b) {
    tPass <- numeric(iterations)
    cramersVPass <- numeric(iterations)
    sampleSizeVec <- numeric(iterations)
    
    for (j in 1:iterations) {
      sampleRaw <- numeric()
      cramersV <- numeric(T)
      
      for (t in 1:T) {
        if (t == 1) {
          sampleRaw <- sample(1:b, n_t[1], replace = TRUE, prob = probabilities)
        } else {
          sampleRaw <- c(sampleRaw, sample(1:b, n_t[t] - n_t[t - 1], replace = TRUE, prob = probabilities))
        }
        
        repVal <- rep(n_t[t] / b, b)
        sample <- table(factor(sampleRaw, levels = 1:b))
        chiSquare <- sum(((sample - repVal)^2) / repVal)
        cramersV[t] <- sqrt(chiSquare / n_t[t]) / sqrt(b - 1)
      }
      
      cramersVUpper <- sqrt(chiSquareU / n_t) / sqrt(b - 1)
      
      if (any(cramersV > cramersVUpper)) {
        tPass[j] <- min(which(cramersV > cramersVUpper))
        cramersVPass[j] <- cramersV[tPass[j]]
        sampleSizeVec[j] <- n_t[tPass[j]]
      } else {
        tPass[j] <- NaN
        cramersVPass[j] <- NaN
        sampleSizeVec[j] <- NaN
      }
    }
    
    list(
      tPass = tPass,
      cramersVPass = cramersVPass,
      sampleSizeVec = sampleSizeVec
    )
  }
  
  # Initialize result vectors
  results <- list(
    pVec = c(), n1Vec = c(), TVec = c(), tPass = c(), cramersVVec = c(),
    powerVec = c(), sampleSize = c(), falsePositive = c()
  )
  counter <- 0
  
  for (m in seq_along(TVecAcceptable)) {
    T <- TVecAcceptable[m]
    n_1 <- n1VecAcceptable[m]
    s <- (1 - n_1) / T
    n_t <- round(((1:T) * n_1) + (((1:T) - 1) * ((1:T) + 2) / 2) * s)
    
    # Run simulations for True Positive
    for (i in 1:simulations) {
      counter <- counter + 1
      p <- runif(1, 0, 1)
      probabilities <- sapply(bins, truncatedGeometricPMF, p = p, b = b)
      
      experiment_result <- simulate_one_experiment(probabilities, T, n_t, chiSquareU, iterations, b)
      
      results$pVec[counter] <- p
      results$n1Vec[counter] <- n_1
      results$TVec[counter] <- T
      results$tPass[counter] <- median(experiment_result$tPass, na.rm = TRUE)
      results$cramersVVec[counter] <- median(experiment_result$cramersVPass, na.rm = TRUE)
      results$sampleSize[counter] <- median(experiment_result$sampleSizeVec, na.rm = TRUE)
      results$powerVec[counter] <- sum(!is.na(experiment_result$tPass)) / iterations
      results$falsePositive[counter] <- NaN
    }
    
    # Run simulations for false positives (null hypothesis true)
    for (i in 1:simulations) {
      counter <- counter + 1
      probabilities <- rep(1 / b, b)
      
      experiment_result <- simulate_one_experiment(probabilities, T, n_t, chiSquareU, iterations, b)
      
      results$pVec[counter] <- NaN
      results$n1Vec[counter] <- n_1
      results$TVec[counter] <- T
      results$tPass[counter] <- NaN
      results$cramersVVec[counter] <- median(experiment_result$cramersVPass, na.rm = TRUE)
      results$sampleSize[counter] <- NaN
      results$powerVec[counter] <- NaN
      results$falsePositive[counter] <- sum(!is.na(experiment_result$tPass)) / iterations
    }
    
    # make it a cleaner number for potential more bin count.
    print(paste("Percent complete:", round(m / length(TVecAcceptable) * 100, 2)))
  }
  
  # Combine results into a data frame
  
  data <- data.frame(
    p = results$pVec, n1 = results$n1Vec, T = results$TVec, tPass = results$tPass,
    cramersV = results$cramersVVec, Power = results$powerVec, SampleSize = results$sampleSize,
    falsePositive = results$falsePositive
  )
  
  return(data)
}
