
#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

source("toolFunctions.R")

#### Generate data
## Common variables
set.seed(1969-08-18)

alpha = -3.5
beta = 0.33/45 # Assuming that sd(dbh) = 45, which is more or less the case with lines 29-30
sigmaProc = 3.2

n = 1e4

## Growth function
growth = function(dbh, alpha, beta)
	return(45*exp(alpha + beta*dbh)) # Assuming that sd(dbh) = 45, which is more or less the case with lines 29-30

## Explanatory variable
dbh_start = rgamma(n, shape = 150^2/2000, rate = 150/2000)
dbh_end = dbh_start + rgamma(n, shape = growth(dbh_start, alpha, beta)^2/sigmaProc, rate = growth(dbh_start, alpha, beta)/sigmaProc)

range(dbh_start)
range(dbh_end)

data = data.table(tree_id = rep(1:n, each = 2), dbh = c(rbind(dbh_start, dbh_end)))
data[seq(1, .N, by = 2), growth := dbh_end - dbh_start]

sd_dbh = data[, sd(dbh)]

#### Model
## Common variables
n_chains = 4
n_data = data[, .N]

## Data
stanData = list(
	n = n_data,
	n_indiv = n,
	sd_dbh = sd_dbh,

	parents_index = seq(1, n_data - 1, by = 2),
	children_index = seq(2, n_data, by = 2),
	dbh = data[, dbh]
)

## Compile
model = cmdstan_model("./verifSigmaProc.stan")

## Run
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 250, chains = n_chains)

results$print()

#### Plots
sigmaProc_array = results$draws("sigmaProc")
lazyPosterior(draws = sigmaProc_array, fun = dgamma, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5)

mean(sigmaProc_array)
mean(rgamma(1e6, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5))
sd_dbh^2*mean(sigmaProc_array)
sd_dbh^2*mean(rgamma(1e6, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5))
sigmaProc

curve(dgamma(x, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5), 0, 0.01)

lazyTrace(draws = sigmaProc_array)

set.seed(2)
y = rgamma(1e6, shape = 3.64^2/1.5, rate = 3.64/1.5)
mean(y)
var(y)

set.seed(2)
y_tilde = rgamma(1e6, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5)

mean(abs(y/sd_dbh^2 - y_tilde))

y = rgamma(1e6, shape = 3.64^2/1.5, rate = 3.64/1.5)
y_tilde = rgamma(1e6, shape = 3.64^2/1.5, rate = sd_dbh^2*3.64/1.5)

hist(y/sd_dbh^2 - y_tilde)
