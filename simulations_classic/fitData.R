
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
data = readRDS("./dummyData_moreVariance.rds")

## Saving options
savingPath = "./"

basename = paste0("indiv=", data[["infos"]]["n_indiv"], "_plot=", data[["infos"]]["n_plot"],
	"_measurements=", data[["infos"]]["n_measurements"], "_deltaT=", data[["infos"]]["delta_t"])

#### Fit data
## Common variables
growth_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "growth[:digit:]")]
temperature_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "avg_temperature[:digit:]")]

n_chains = 3

## Stan data
stanData = list(
	# Number of data
	n_indiv = data[["infos"]]["n_indiv"], # Total number of individuals
	n_obs_growth_per_indiv = data[["infos"]]["n_measurements"] - 1, # Total number of growth observations per individual
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
model = cmdstan_model("./growth_classic.stan")

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
results$save_object(file = paste0(basename, "_results_moreVariance.rds"))

####! Crash test zone
# #### Check correlation growth
# cor(data[["treeData"]][, dbh2 - dbh1], data[["treeData"]][, dbh3 - dbh2])
# cor(data[["treeData"]][, dbh2 - dbh1], data[["treeData"]][, dbh4 - dbh3])
# cor(data[["treeData"]][, dbh2 - dbh1], data[["treeData"]][, dbh5 - dbh4])
# cor(data[["treeData"]][, dbh2 - dbh1], data[["treeData"]][, dbh6 - dbh5])

# cor(data[["treeData"]][, dbh3 - dbh2], data[["treeData"]][, dbh4 - dbh3])
# cor(data[["treeData"]][, dbh3 - dbh2], data[["treeData"]][, dbh5 - dbh4])
# cor(data[["treeData"]][, dbh3 - dbh2], data[["treeData"]][, dbh6 - dbh5])

# cor(data[["treeData"]][, dbh4 - dbh3], data[["treeData"]][, dbh5 - dbh4])
# cor(data[["treeData"]][, dbh4 - dbh3], data[["treeData"]][, dbh6 - dbh5])

# cor(data[["treeData"]][, dbh5 - dbh4], data[["treeData"]][, dbh6 - dbh5])

# #### Check if a sum of lognormal random variables can be approximated by one lognormal random variable
# # The package lognorm is based on http://www.m-hikari.com/ams/ams-2013/ams-125-128-2013/39511.html
# # WKB Approximation for the Sum of Two Correlated Lognormal Random Variables
# library(lognorm)

# # generate nSample values of two lognormal random variables
# mu1 = log(110)
# mu2 = log(100)
# mu3 = log(77)
# mu4 = log(29)
# sigma1 = 0.25
# sigma2 = 0.15
# sigma3 = 0.54
# sigma4 = 0.07

# (coefSum = estimateSumLognormal(c(mu1, mu2, mu3, mu4), c(sigma1, sigma2, sigma3, sigma4)))

# X1 = rlnorm(1e6, meanlog = mu1, sdlog = sigma1)
# X2 = rlnorm(1e6, meanlog = mu2, sdlog = sigma2)
# X3 = rlnorm(1e6, meanlog = mu3, sdlog = sigma3)
# X4 = rlnorm(1e6, meanlog = mu4, sdlog = sigma4)

# Y = X1 + X2 + X3 + X4
# plot(density(Y))
# curve(dlnorm(x, meanlog = coefSum["mu"], sdlog = coefSum["sigma"]),
# 	to = 1200, lwd = 2, col = "#A128AF", add = TRUE)

# dens = density(Y)
# errAbs = abs(dens$y - dlnorm(x = dens$x, meanlog = coefSum["mu"], sdlog = coefSum["sigma"]))
# errRel = abs((dens$y - dlnorm(x = dens$x, meanlog = coefSum["mu"], sdlog = coefSum["sigma"]))/dlnorm(x = dens$x, meanlog = coefSum["mu"], sdlog = coefSum["sigma"]))
# ind = which(dens$x < 600) 
# plot(dens$x[ind], errRel[ind])
# plot(dens$x[ind], 100*errRel[ind])
# plot(dens$x[ind], errAbs[ind])
