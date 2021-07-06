
#### Aim of prog: estimate the parameters that generated the data in createDummyData.R

#### Clear memory and load packages
rm(list = ls())

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### Create data
source("createDummyData.R")

#### Stan model
## Compile
# model = stan_model(file = "./hiddenMarkov.stan")

## Define stan variables
maxIter = 2e3

stanData = list(N = data[, .N] - n,
	K = 2,
	Y = data[!is.na(incr_m), incr_m])

# results = sampling(object = model, data = stanData,
# 	iter = maxIter,
# 	control = list(adapt_delta = 0.85, max_treedepth = 11),
# 	include = FALSE)

results = stan(file = "hiddenMarkov.stan", data = stanData,
	iter = maxIter)

print(results, pars = "lambda")
print(results, pars = "gamma")

traceplot(results)
