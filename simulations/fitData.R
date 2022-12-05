
#### Aim of prog: Fit the dummy data generated from the script "createData.R"

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
source("./toolFunctions.R")

#### Load data
## General data (list)
data = readRDS("./dummyData.rds")

## Saving options
savingPath = "./"

basename = paste0("indiv=", data[["infos"]]["n_indiv"], "_plot=", data[["infos"]]["n_plot"],
	"_measurements=", data[["infos"]]["n_measurements"], "_deltaT=", data[["infos"]]["delta_t"])

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
	dbh_init = unlist(data[["treeData"]][, "obs_dbh1"]),

	# Explanatory variables
	sd_dbh = data[["scaling"]]["sd_dbh_orig"], # To standardise the initial dbh

	tas = data[["treeData"]][, ..temperature_cols], # Temperature
	tas_mu = data[["scaling"]]["mu_temp"], # To centre the temperature
	tas_sd = data[["scaling"]]["sd_temp"], # To standardise the temperature

	# Provided parameters
	beta0 = data[["parameters"]]["beta0"],
	beta1 = data[["parameters"]]["beta1"],
	beta2 = data[["parameters"]]["beta2"],

	beta3 = data[["parameters"]]["beta3"],
	beta4 = data[["parameters"]]["beta4"],

	# Observation error (not yet simulated in the data actually, just used for the likelihood)
	sigmaObs = data[["sigmaObs"]]
)

## Compile model
model = cmdstan_model("./growth.stan")

## Run model
start_time = Sys.time()

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 750, iter_sampling = 1500, save_warmup = TRUE,
	max_treedepth = 12, adapt_delta = 0.9)

end_time = Sys.time()

results$cmdstan_diagnose()

#### Plots
## Lazy trace plots
lazyTrace(draws = results$draws("beta0"), val1 = data[["parameters"]]["beta0"], main = paste("beta0 =", data[["parameters"]]["beta0"]))
lazyTrace(draws = results$draws("beta1"), val1 = data[["parameters"]]["beta1"], main = paste("beta1 =", data[["parameters"]]["beta1"]))
lazyTrace(draws = results$draws("beta2"), val1 = data[["parameters"]]["beta2"], main = paste("beta2 =", data[["parameters"]]["beta2"]))
lazyTrace(draws = results$draws("beta3"), val1 = data[["parameters"]]["beta3"], main = paste("beta3 =", data[["parameters"]]["beta3"]))
lazyTrace(draws = results$draws("beta4"), val1 = data[["parameters"]]["beta4"], main = paste("beta4 =", data[["parameters"]]["beta4"]))
lazyTrace(draws = results$draws("sigmaProc"), val1 = data[["parameters"]]["sigmaProc"], main = "sigmaProc")
lazyTrace(draws = results$draws("latent_dbh_parents[1]"),
	val1 = data[["treeData"]][1, dbh1]/data[["scaling"]]["sd_dbh_orig"])

results$save_output_files(dir = savingPath, basename = paste0(basename, "_diagnose"), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(basename, "_results.rds"))

indiv = 1

for (i in 1:(data[["infos"]]["delta_t"] - 1))
{
	current_dbh = paste0("dbh", i)
	next_dbh = paste0("dbh", i + 1)
	lazyPosterior(draws = data[["scaling"]]["sd_dbh_orig"]*results$draws(paste0("latent_growth[", indiv, ",", i, "]")), fun = NULL,
		val1 = data[["treeData"]][indiv, ..next_dbh] - data[["treeData"]][indiv, ..current_dbh])
}
