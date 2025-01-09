# AMEquivalency

R package used to determine equivalency between reference and candidate process. A big thank you to Cesar Gomez Guillen (https://github.com/cfgomezguil) for helping to optimize this code. I also want to thank Brenda Valdez Mesta for testing the code and for providing helpful insights.    

## How to install

Install and load the package devtools and then type the following into your R script:

```
install_github("colinmichaellynch/AMEquivalency")
library(AMEquivalency)
```

## Tutorial 

The R script tutorial.R was created to demonstrate how the AMEquivalency package is meant to be utilized. This script walks through the steps highlighted in the manuscript "Determining univariate equivalency of additively
manufactured parts" (Lynch et al., in prep). It first characterizes a reference distribution (ReferenceData.csv) and tests its stability. It then uses simulations to determine what the optimal sequential sampling strategy is for determining whether a distribution created by some candidate process is equivalent to the distribution produced by the reference process. Finally, it provides tests whether a sample is equivalent to the reference. If it is equivalent, then the experimentor is expected to continue sampling. If they are not equivalent, then the experimentor stops sampling. There are an additional two example datasets available which act as candidate distributions that is currently undergoing sequential sampling. One of these examples (CandidateEquivalentSample.csv) is equivalent, and so the test will tell the user to keep sampling. The other example (CandidateNotEquivalentSample.csv) is not equivalent, so the test will tell the user to stop sampling. 

You must download tutorial.R, ReferenceData.csv, CandidateEquivalentSample.csv, and CandidateNotEquivalentSample.csv independently of the package. Make sure that all files are in the same working directory. These files are available in the tutorial folder. 

Note that data files need to have a specific structure to be process properly. The reference data will have n rows and m columns, where each column represents n samples from each subgroup:

![referenceData](https://github.com/user-attachments/assets/8d721a17-b171-4836-b30b-6827321d17af)

Conversely, the candidate dataset only needs a single column with the values of all observations: 
![candidateData](https://github.com/user-attachments/assets/a4de25a1-46cb-4dc2-9603-aed92818a7c0)
