
#### Clear memory and load packages
rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### Generate dummy data
## Common variables
n = 200 # Number of individuals

set.seed(1969-08-18) # Woodstock seed

latentMean = 5
latentVar = 12

latentShape = latentMean^2/latentVar
latentRate = latentMean/latentVar

latent = rgamma(n, shape = latentShape, rate = latentRate)

(mean(latent))
(var(latent))

measurementError = 2

measuredDbh = rnorm(n, mean = latent, sd = measurementError)

#### Stan model
## Compile
# model = stan_model(file = "./hiddenMarkov.stan")

## Define stan variables
maxIter = 1e4

stanData = list(N = n,
	Y = measuredDbh)

results = stan(file = "latent_transformed.stan", data = stanData,
	iter = maxIter, control = list(adapt_delta = 0.99))

print(results, pars = "latentMean")
print(results, pars = "latentVar")
print(results, pars = "measurementError")

traceplot(results, pars = c("latentShape", "latentRate", "measurementError"))
traceplot(results, pars = c("latentMean", "latentVar", "measurementError"))

## Comparison true dbh estimation versus latent
# trueDbh[190]        8.33
# trueDbh[191]        3.12
# trueDbh[192]       17.37
# trueDbh[193]        4.19
# trueDbh[194]        1.24
latent[190:194]

# trueDbh[195]        1.47
# trueDbh[196]        1.63
# trueDbh[197]        9.68
# trueDbh[198]        6.27
# trueDbh[199]        2.20
# trueDbh[200]        5.50
tail(latent)

estimatedTrueDBH = get_posterior_mean(results, pars = paste0("trueDbh[",1:200,"]"))[, "mean-all chains"]

var(latent - estimatedTrueDBH)
hist(latent - estimatedTrueDBH)
