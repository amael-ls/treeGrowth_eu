
# This prog generates data for toy.stan, toyPara_GPUs.stan and toyPara_reduce_sum.stan

####
rm(list = ls())
graphics.off()

library(data.table)

options(max.print = 500)

#### Variables
# Number of individuals
n_indiv = 200

# Set seed
set.seed(1969-08-18) # Woodstock seed

# Number of states per individual
nbYearsPerIndiv = sample(x = 4:8, size = n_indiv, replace = TRUE)
barplot(table(nbYearsPerIndiv))
