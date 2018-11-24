# BWMR
BWMR (Bayesian Weighted Mendelian Randomization), is an efficient statistical method to infer the causality between a risk exposure factor and a trait or disease outcome, based on GWAS summary statistics. 'BWMR' package provides the estimate of causal effect with its standard error and the P-value under the test of causality.


# Installation
To install the development version of BWMR, it's easiest to use the 'devtools' package.
'''
# install.packages("devtools")
library(devtools)
install_github("jiazhao97/BWMR")
'''


# Usage
The main function is *BWMR*. You can find an example by running
'''
library(ggplot2)
library(BWMR)
example(BWMR)
'''


# Development
This R package is developed by Jia Zhao.

