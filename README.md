## Circadian Analysis: Rhythmicity detection and DE analysis 

Readme.txt file for Rhythmicity and DE Analysis created 1/24/2019

# Required R packages: 
minpack.lm v1.2-1
Snowfall v1.84-6.1
doParallel v1.014
Lme4 v1.1-19
Value v2.14.1
This code has been tested on R.3.5.1 

# Installation Instructions: 
Download and unzip the folders to the desktop 
"rhythmicity" folder contains code for rhythmicity detection
"DE" folder contains code for DE analysis
Save each folder to the desktop
Install required packages: 
Example:
```
install.packages("minpack.lm") 
library(minpack.lm)
```
Installation time: < 2 min

# Demo:
A subset of expression, and clinical data for 104 control subjects is available in csv form in the rhythmicity folder. 
Clinical variables of interest include:
corrected time of death (TOD) in ZT. 
Hospital site (Mount Sinai vs Pitt)
Diagnosis (control vs schizophrenia)
Gender 

# Expression data includes: 
13914 genes (all genes that were used in the full analysis in the paper)
Data has already been filtered according to method discussed in the supplement 
Expression units are in log2 CPM 

fitSinCurve.R and Curve_Drawing.R are external functions required to run the main analysis in RhythmicityCode.R 
BestModelSelection.R is an external function required to run the main DE.R analysis
Run all external functions, so that they are in the R environment 
Run code as is in RhythmicityCode.R 
Run code as is in DE.R 
*NOTE: Each method requires permutation. Permutations were set to 10 to save time and computation space. Permuted files are included in the DE and Rhythmicity folders to save reviewers time. If permutation files are recreated by the reviewer, output will be slightly different due to the randomness of permutation.

# Expected output:
Example_result.csv is the example output from the observed_para_c_sorted variable. 
PDF plots are available in the Rhythmicity folder 
Example_result2.csv is the example output for the shift in rhythmicity analysis 
*NOTE: Because only control subjects are used in the example data, two groups were created using age to demonstrate gain and loss of rhythmicity 

Example_result3.csv in the DE folder includes 
*NOTE: Because only control subjects are used in the example data, effect of a binary indicator for time of day (morning vs night) is analyzed. 

Estimated run time: 
30min-1hr

