# BWMR
BWMR (Bayesian Weighted Mendelian Randomization), is an efficient statistical method to infer the causality between a risk exposure factor and a trait or disease outcome, based on GWAS summary statistics. 'BWMR' package provides the estimate of causal effect with its standard error and the P-value under the test of causality.


# Installation
To install the development version of BWMR, it's easiest to use the 'devtools' package.
```
# install.packages("devtools")
library(devtools)
install_github("jiazhao97/BWMR")
```


# Usage
[The 'BWMR' vignette](https://github.com/jiazhao97/BWMR/blob/master/vignettes/BWMR_package.pdf?raw=true) will provide a good start point for Mendelian randomization analysis using BWMR package. The main function is BWMR. You can find examples by running
```
library(ggplot2)
library(BWMR)
example(BWMR)
```


# Reference
Jia Zhao, Jingsi Ming, Xianghong Hu, Gang Chen, Jin Liu, Can Yang, Bayesian weighted Mendelian randomization for causal inference based on summary statistics, Bioinformatics, btz749, https://doi.org/10.1093/bioinformatics/btz749


# Development
This R package is developed by Jia Zhao.

