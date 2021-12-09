
#### Aim of prog: Fit the dummy data created with createData.R with a hierarchical model (in stan). Check the residuals

#### Clear memory and load libraries
rm(list = ls())
graphics.off()

library(data.table)
library(cmdstanr)
library(viridis)
library(DHARMa)

#### Data
## Load data
dt = readRDS("./data.rds")

## Create stan data
n_data = dt[, .N]
n_group = length(unique(dt[, group]))
data_per_group = dt[, .N, by = group][, N]

stanData = list(
	n_data = n_data,
	n_group = n_group,
	data_per_group = data_per_group,
	x = dt[, x],
	y = dt[, y]
)

#### Run
## Parameters
maxIter = 2e3
n_chains = 3

## Compile model
model = cmdstan_model("./model.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 10, adapt_delta = 0.8)

diagnose = results$cmdstan_diagnose()

#### Check the fit
estimatedParameters = c(beta0 = mean(results$draws("beta0")),
	sigma_beta = mean(results$draws("sigma_beta")),
	beta0_1 = 0,
	beta0_2 = 0,
	beta0_3 = 0,
	beta0_4 = 0,
	beta1 = mean(results$draws("beta1")),
	sigma_res = mean(results$draws("sigma_res")))
for (i in 1:n_group)
	estimatedParameters[paste0("beta0_", i)] = mean(results$draws(paste0("beta0_group[", i, "]")))

realParameters = readRDS("./realParameters.rds")
print(realParameters)
print(estimatedParameters)

#### Analysing residuals
