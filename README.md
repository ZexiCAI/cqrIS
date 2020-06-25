
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cqrIS: Censored Quantile Regression with Induced Smoothing

<!-- badges: start -->

<!-- badges: end -->

The goal of cqrIS is to fit a censored quantile regression model to
ordinary survival data. It also provides tools to fit a censored
quantile regression model to functional covariates and a recurrent event
model to doubly-censored recurrent events. Improvements on interquantile
smoothness are also available in this package.

## Installation

You can install the released version of cqrIS from GitHub with:

``` r
devtools::install_github("ZexiCAI/cqrIS")
library(cqrIS)
```

## Example

This is a basic example which shows you how to fit a censored quantile
regression model to ordinary survival data, and obtain its smoothed
version:

``` r
library(cqrIS)
dat <- ordin.sam200.cen25.homo
res.PH <- estPH(Z=dat[,-c(1:2)], X=dat[,1], cen=dat[,2])
est.PH <- res.PH[[1]]

res.cqrIS <- cqrIS(Z=dat[,-c(1:2)], X=dat[,1], cen=dat[,2])
est.cqrIS <- res.cqrIS[[1]]
```

See the package vignette for more information and detailed instruction.

## Reference

  - Cai, Z. and Sit, T. (2020+), “Censored Quantile regression with
    Induced Smoothing,” Working Paper.

  - Jiang, F., Cheng, Q., Yin, G. and Shen, H. (2020), “Generalizing
    Quantile Regression for Counting Processes with Applications to
    Recurrent Events,” *Journal of the American Statistical
    Association*, **115**, 931-944.

  - Peng, L. and Huang, Y. (2008), “Survival Analysis with Quantile
    Regression Models,” *Journal of the American Statistical
    Association*, **103**, 637-649.

  - Sun, X., Peng, L., Huang, Y. and Lai, H.J. (2016), “Generalizing
    Quantile Regression for Counting Processes with Applications to
    Recurrent Events,” *Journal of the American Statistical
    Association*, **111**, 145-156.
