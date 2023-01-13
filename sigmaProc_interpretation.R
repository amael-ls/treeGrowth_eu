
#### Aim of prog: Comparing sigmaProc (SSM approach) vs sigmaProc (classic approach)
# Comments:
# The package lognorm is based on http://www.m-hikari.com/ams/ams-2013/ams-125-128-2013/39511.html
#	"WKB Approximation for the Sum of Two Correlated Lognormal Random Variables".
#	I use this package to approximate the sum of lognormal random variables by a single lognormal distribution. The article looks fine,
#	although the author did not check the relative error of his method, only the absolute error. The relative error shows that his method
#	is unable to correclty fit the tail or the head of the distribution. Other methods have been developed where the precision in the tail
#	or head is also discussed in https://ieeexplore.ieee.org/document/1578407 "A Flexible Lognormal Sum Approximation Method".
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(lognorm)
