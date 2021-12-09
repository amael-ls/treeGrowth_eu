
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
n_group = length(unique(dt[, group]))

stanData = list()

#### Run
## Parameters
maxIter = 2e3
n_chains = 3

## Compile model
model = cmdstan_model("./model.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 10, adapt_delta = 0.8)

## Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-", time_ended, ".rds"))

diagnose = results$cmdstan_diagnose()

#### Analysing residuals
