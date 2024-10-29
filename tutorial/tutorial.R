#In this tutorial, we will walk through the steps necessary to first characterize a reference distribution and then compare this distribution to some candidate distribution to determine whether or not the two are equivalent. It will cover sequential sampling, the required data format, and other topics related to equivalency. We recommend that you first read "Determining univariate equivalency of additively manufactured parts" by Lynch et al., in prep to understand the broader framework for this methodology.

#if you have not yet done so already, you will need to clear the workspace and then install AMEquivalency, which will require a function from the library devtools:

rm(list=ls())

library(devtools)
install_github("colinmichaellynch/AMEquivalency")
library(AMEquivalency)

#Next, we will set our experimental constraints: 

b = 8 #number of bins
alpha = .01 #significance level used in stopping rule
n_max = 499 #highest possible sample size acceptable in study
simulations = 100 #number values of rho tested in simulations to find optimal sampling strategy
iterations = 100 #number of simulations per parameter value of rho 
N_1T_Max = 40 #highest possible initial sample size and number of sequential samples
wVec = c(1/4, 1/4, 1/2) #weight vector for design optimization. The first element corresponds to power, the second to sample size, the last to the false positive rate
L1 = 0 #lowest acceptable power level
T1 = 1 #highest acceptable power level
T2 = 8 #lowest acceptable sample size
U2 = 41 #highest acceptable sample size
T3 = 0 #lowest acceptable false positive rate
U3 = 1 #highest acceptable sample size
k = 2 #concavity of desirability functions 

#Now we will characterize the reference dataset. Here, we will play around with an artificially created dataset for illistrative purposes. Feel free to replace this dataset with your own. This dataset must be n by m in size, where n is the samples per subgroups and m is the number of subgroups. Each cell must contain a numeric value. The minimal sample size is n = 2 and m = 20. 

dataRef = read.csv("ReferenceData.csv", header = FALSE)

#we will first determine whether or not the process is stable using an x-bar chart. This function will also count the number of points that fall outside of the control limits
outOfControlCount = processStability(dataRef, show.plot = TRUE)
print(paste("# Out of Control Points:", outOfControlCount))

#If the process is stable, then we can characterize the distribution by dividing it into b equally-likely bins with using percentiles. Ensure that the first and last bin limits are set to negative infinity and infinity, respectively  
dataRef = as.vector(as.matrix(dataRef))
percentileVec = seq(0, 1, length.out = b+1)
bins = quantile(dataRef, percentileVec)
bins[1] = -Inf
bins[b+1] = Inf

#We will use sequential sampling to measure the candidate process. We first set the maximal initial sample size/number of samples we are willing to take. This then gives us a family of possible strategies we can use to sample with. 

acceptableStrats = possibleSamplingStrats(N_1T_Max)
TVecAcceptable = acceptableStrats[[1]]
n1VecAcceptable = acceptableStrats[[2]]

#given these possible sampling strategies (each pair of elements in TVecAcceptable and n1VecAcceptable represent a possible strategy), we can run simulations to determine the 3 performance metrics of our designs: power of each design, the false positive rate, and the expected sample size given that the two processes are not equivalent. The following function runs simulations across many situations where the two processes are equivalent (rho = 0), are slightly not equivalent (small rho), or are very much not equivalent (big rho) 

dataSimulated = simulateExperiment(b, alpha, n_max, simulations, iterations, TVecAcceptable, n1VecAcceptable)

#We can then use desirability functions to determine which strategy best trades off our 3 performance metrics 

dataOptimal = optimalStrategy(dataSimulated, wVec, L1, T1, T2, U2, T3, U3, k)

print(paste("Power =", dataOptimal[[1]], ": Sample Size (non-equivalent) =", dataOptimal[[2]], ": False Positive Rate =", dataOptimal[[3]]))

N_1 = dataOptimal[[4]]
T = dataOptimal[[5]]

#Now that we have our optimal design, we use this function to figure out what the sequence of sample sizes we need to take. This is in terms of the cumulative sample size. For example, if the first sample is of size 5, and the second of size 9, that just means that you need to include an additional 4 samples before you can run another test

sequence = cumulativeSamplingSequence(N_1, T)
print(sequence)

#Now we will do an equivalency test on a sample which is the equivalent to the reference. In our example case, we are assuming that we are already partway through our sequential sampling regimen. We will first import the candidate data: 

dataCandEquivalent = read.csv("CandidateEquivalentSample.csv", header = FALSE)

#Then we will measure its chi square value (how different the sample is from the expectation that all observations should be uniformly distributed among the reference bins). The graph illustrates how well this uniform assumption holds

chiSquare = candidateChiSquared(dataCandEquivalent, bins, show.plot = TRUE)

#we then calculate cramers v (which normalizes the chi square given the current sample size)
N = length(dataCandEquivalent)
cramersV = candidateCramersV(chiSquare, b, N)

#and then we test whether or not the cramer's v is within its expected bounds given that the processes are equivalent. In this first case, the two processes are equivalent, so we get a 'false' value from this function: 
stoppingRule(cramersV, N, b, alpha)

#in this case, you would keep sampling until you've reached the end of the sequence as determined by cumulativeSamplingSequence (assuming that every sample says that the two processes are equivalent)

#we now go through the same process for the process which is not equivalent to the reference. Here, the final decision value is 'true'. If this happens, you can stop and determine that the two are not equivalent. 
dataCandNotEquivalent = read.csv("CandidateNotEquivalentSample.csv", header = FALSE)

chiSquare = candidateChiSquared(dataCandNotEquivalent, bins, show.plot = TRUE)
N = length(dataCandNotEquivalent)
cramersV = candidateCramersV(chiSquare, b, N)
stoppingRule(cramersV, N, b, alpha)

