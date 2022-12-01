
#### Aim of prog: Fit the dummy data generated from the script "createData.R"

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Load data
## General data (list)
data = readRDS("./dummyData.rds")

#### Fit data
## Common variables
growth_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "growth[:digit:]")]
temperature_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "temperature[:digit:]")]

n_chains = 3

## Stan data
stanData = list(
	# Number of data
	n_indiv = data[["infos"]]["n_indiv"], # Total number of individuals
	n_obs_growth = data[["infos"]]["n_measurements"] - 1, # Total number of growth observations
	n_latent_growth = data[["infos"]]["n_annual_growth_per_indiv"], # Number of growing years (for latent state)
	delta_t = data[["infos"]]["delta_t"], # Number of growing years (for latent state)

	# Observations
	avg_yearly_growth_obs = data[["treeData"]][, ..growth_cols],

	# Explanatory variables
	sd_dbh = data[["scaling"]]["sd_dbh_orig"], # To standardise the initial dbh

	tas = data[["treeData"]][, ..temperature_cols], # Temperature
	tas_mu = data[["scaling"]]["mu_temp"], # To centre the temperature
	tas_sd = data[["scaling"]]["sd_temp"] # To standardise the temperature
)

## Compile model
model = cmdstan_model("./growth.stan")

## Run model
start_time = Sys.time()

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 1500, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.9)

end_time = Sys.time()
