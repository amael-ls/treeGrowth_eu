
#### Aim of prog: To sum-up the pre-runs, such as mean and sd of estimated parameters among species and within species (if run > 1)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

