
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
delta_t = 12 # Choose in accordance with the generated data sets from createData
n_indiv = 400
n_plot = 50

data = readRDS(paste0("dummyData_plot=", n_plot, "_indiv=", n_indiv, "_deltaT=", delta_t, ".rds"))

## Saving options
savingPath = "./"

basename = paste0("indiv=", data[["infos"]]["n_indiv"], "_plot=", data[["infos"]]["n_plot"], "_deltaT=", data[["infos"]]["delta_t"])

#### Fit data
## Common variables
growth_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "growth[:digit:]")]
temperature_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "temperature[:digit:]")]
temperature_cols = temperature_cols[!stri_detect(temperature_cols, regex = "avg_temperature[:digit:]")]

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
model = cmdstan_model("./growth.stan", cpp_options = list(stan_opencl = TRUE))

## Run model
start_time = Sys.time()

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 750, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 10, adapt_delta = 0.8, opencl_ids = c(0, 0))

end_time = Sys.time()

print(paste0("Running time: ", round(end_time - start_time, 2)))

results$cmdstan_diagnose()

#### Plots
## Lazy trace plots
# lazyTrace(draws = results$draws("beta0"), val1 = data[["parameters"]]["beta0"], main = paste("beta0 =", data[["parameters"]]["beta0"]))
# lazyTrace(draws = results$draws("beta1"), val1 = data[["parameters"]]["beta1"], main = paste("beta1 =", data[["parameters"]]["beta1"]))
# lazyTrace(draws = results$draws("beta2"), val1 = data[["parameters"]]["beta2"], main = paste("beta2 =", data[["parameters"]]["beta2"]))
# lazyTrace(draws = results$draws("beta3"), val1 = data[["parameters"]]["beta3"], main = paste("beta3 =", data[["parameters"]]["beta3"]))
# lazyTrace(draws = results$draws("beta4"), val1 = data[["parameters"]]["beta4"], main = paste("beta4 =", data[["parameters"]]["beta4"]))
# lazyTrace(draws = results$draws("sigmaProc"), val1 = data[["parameters"]]["sigmaProc"], main = "sigmaProc")
# lazyTrace(draws = results$draws("latent_dbh_parents[1]"),
# 	val1 = data[["treeData"]][1, dbh1]/data[["scaling"]]["sd_dbh_orig"],
# 	val2 = data[["treeData"]][1, obs_dbh1]/data[["scaling"]]["sd_dbh_orig"],
# 	label1 = "real",
# 	label2 = "obs")

basename = paste0("sigmaObs=0.1_", basename)
results$save_output_files(dir = savingPath, basename = paste0(basename, "_diagnose"), timestamp = FALSE, random = FALSE)
results$save_object(file = paste0(basename, "_results.rds"))

# indiv = 1

# for (i in seq_len(data[["infos"]]["n_annual_growth_per_indiv"] - 1))
# {
# 	current_dbh = paste0("dbh", i)
# 	next_dbh = paste0("dbh", i + 1)
# 	lazyPosterior(draws = data[["scaling"]]["sd_dbh_orig"]*results$draws(paste0("latent_growth[", indiv, ",", i, "]")), fun = NULL,
# 		val1 = data[["treeData"]][indiv, ..next_dbh] - data[["treeData"]][indiv, ..current_dbh])
# }

# indiv = 100

# for (i in seq_along(data[["infos"]]["delta_t"] - 1))
# {
# 	current_dbh = paste0("dbh", i)
# 	next_dbh = paste0("dbh", i + 1)
# 	lazyPosterior(draws = data[["scaling"]]["sd_dbh_orig"]*results$draws(paste0("latent_growth[", indiv, ",", i, "]")), fun = NULL,
# 		val1 = data[["treeData"]][indiv, ..next_dbh] - data[["treeData"]][indiv, ..current_dbh])
# }
