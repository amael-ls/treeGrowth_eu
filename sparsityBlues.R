
## This script is to follow the blog post from Michael Bettancourt, that can be found at:
# https://betanalpha.github.io/assets/case_studies/modeling_sparsity.html#1_Fading_into_Irrelevance

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Generate data
