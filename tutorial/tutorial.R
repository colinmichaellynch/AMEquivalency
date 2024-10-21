### --- prepare workspace --- ### 

#clear directory, load packages, load functions
#https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

rm(list=ls())

setwd("/path/to/directory")

source("processStability.R")
source("cumulativeSamplingSequence.R")
source("candidateChiSquared.R")
source("candidateCramersV.R")
source("stoppingRule.R")
source("possibleSamplingStrats.R")
source("truncatedGeometricPMF.R")
source("simulateExperiment.R")
source("optimalStrategy.R")

library(tidyverse)

### set constants

b = 8 #number of bins
alpha = .01 #significance level used in stopping rule
n_max = 499 #highest possible sample size acceptable in study
simulations = 10 #number values of rho tested in simulations to find optimal sampling strategy
iterations = 10 #number of simulations per parameter value rho 
N_1T_Max = 40 #highest 

dataRef = read.csv("exampleData.csv")
dataRef = dataRef[, -1]

save(exampleReferenceData, file="exampleReferenceData.RData")

#check if there are any out of control points
outOfControlCount = processStability(dataRef, show.plot = TRUE)
print(paste("# Out of Control Points:", outOfControlCount))

#optimal sampling strategy

acceptableStrats = possibleSamplingStrats(N_1T_Max)
TVecAcceptable = acceptableStrats[[1]]
n1VecAcceptable = acceptableStrats[[2]]

dataSimulated = simulateExperiment(b, alpha, n_max, simulations, iterations, TVecAcceptable, n1VecAcceptable)

wVec = c(1/4, 1/4, 1/2)
L1 = 0
T1 = 1
T2 = 8 
U2 = 41 
T3 = 0 
U3 = 1 
k = 2

dataOptimal = optimalStrategy(dataSimulated, wVec, L1, T1, T2, U2, T3, U3, k)

print(paste("Power =", dataOptimal[[1]], ": Sample Size (non-equivalent) =", dataOptimal[[2]], ": False Positive Rate =", dataOptimal[[3]]))

N_1 = dataOptimal[[4]]
T = dataOptimal[[5]]

#final sampling strategy 
sequence = cumulativeSamplingSequence(N_1, T)

#test for equivalency 

dataRef = as.vector(as.matrix(dataRef))
percentileVec = seq(0, 1, length.out = b+1)
bins = quantile(dataRef, percentileVec)

dataCandEquivalent = read.csv("CandidateEquivalentSample.csv")
dataCandEquivalent = dataCandEquivalent[, -1]

chiSquare = candidateChiSquared(dataCandEquivalent, bins, show.plot = TRUE)
N = length(dataCandEquivalent)
cramersV = candidateCramersV(chiSquare, b, N)
stoppingRule(cramersV, N, b, alpha)

dataCandNotEquivalent = read.csv("CandidateNotEquivalentSample.csv")
dataCandNotEquivalent = dataCandNotEquivalent[, -1]

chiSquare = candidateChiSquared(dataCandNotEquivalent, bins, show.plot = TRUE)
N = length(dataCandNotEquivalent)
cramersV = candidateCramersV(chiSquare, b, N)
stoppingRule(cramersV, N, b, alpha)
