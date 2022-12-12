
#### Aim of prog: Replicate the small simulation done in Aubry-Kientz et. al (2017) to verify if Eitzel 2013 can be applied for Δt >= 5

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

source("toolFunctions.R")
#### Tool function
## Initiate Y_gen with reasonable values (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("n_indiv", "n_years_growth", "annual_dbh")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide n_indiv, n_years_growth, annual_dbh")

	n_indiv = providedArgs[["n_indiv"]]
	n_years_growth = providedArgs[["n_years_growth"]]
	annual_dbh = providedArgs[["annual_dbh"]]

	latent_growth_gen = matrix(data = NA, nrow = n_indiv, ncol = n_years_growth)

	for (j in 1:n_years_growth)
	{
		current_dbh = paste0("dbh", j)
		next_dbh = paste0("dbh", j + 1)
		annual_growth = annual_dbh[[next_dbh]] - annual_dbh[[current_dbh]]
		latent_growth_gen[, j] = rgamma(n_indiv, annual_growth^2/0.5, annual_growth/0.5) # Average = annual_growth, variance = 0.5
	}

	return(list(latent_yearly_growth = latent_growth_gen))
}

#### Simulate data
## Common variables
set.seed(1969-08-18) # Woodstock seed

n = 1000
n_years = 25

init_dbh = rnorm(n = n, mean = 50, sd = 10)
if (any(init_dbh < 0))
	warning("There are some negative starting dbh")

temperature = rnorm(n = n_years, mean = 7, sd = 2.5)

## Parameters (that should be recovered via Stan)
beta0 = 2
beta1 = 0.02
beta2 = 1.5

sigmaProc = 1.5

## Data
treeData = data.table(tree_id = 1:n, dbh1 = init_dbh)
treeData[, paste0("dbh", 2:n_years) := 0]

for (i in 2:n_years)
{
	current_dbh = paste0("dbh", i)
	previous_dbh = paste0("dbh", i - 1)

	noise = rnorm(n, mean = 0, sd = sigmaProc)
	treeData[, c(current_dbh) := (1 + beta1)*treeData[[previous_dbh]] + beta0 + beta2*temperature[i - 1] + noise]
}

#### Stan model
## Compile model
model = cmdstan_model("Aubry-Kientz2017.stan")

## Stan data
# Common variables
freq_obs = 2
n_chains = 3

# Subset data
keptCols = c("tree_id", paste0("dbh", seq(1, n_years, by = freq_obs)))
treeData_obs = treeData[, ..keptCols]

n_obs_year = length(keptCols) - 1 # -1 because I kept tree_id, which is not an observation

# Compute the growth observations (averaged annual growth)
counter = 1
for (j in 3:length(keptCols))
{
	current_dbh = keptCols[j]
	previous_dbh = keptCols[j - 1]
	growthCol = paste0("growth", counter)

	treeData_obs[, c(growthCol) := (treeData_obs[[current_dbh]] - treeData_obs[[previous_dbh]])/freq_obs]
	counter = counter + freq_obs
}
growthCols = colnames(treeData_obs)[stri_detect(colnames(treeData_obs), regex = "growth[:digit:]")]

# Stan list
stanData = list(
	# Number of data
	n_indiv = n,
	n_years = n_years,
	n_obs_years = n_obs_year,

	# Observations
	dbh_init = treeData_obs[, dbh1],
	observed_averaged_annual_growth = treeData_obs[, ..growthCols],
	freq_obs = freq_obs,

	# Predictors
	temperature = temperature
)

## Initialisation latent states
initVal_Y_gen = lapply(1:n_chains, init_fun, n_indiv = n, n_years_growth = n_years - 1, annual_dbh = treeData)

## Run model
model = cmdstan_model("Aubry-Kientz2017.stan")
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 500, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

results$save_output_files(dir = "./", basename = "Aubry-Kientz2017", timestamp = FALSE, random = TRUE)
results$save_object(file = "./Aubry-Kientz2017_results.rds")

results$print(c("lp__", "beta0", "beta1", "beta2"), max_rows = 20)
results$print(c("lp__", "beta0", "beta1", "beta2", "sigmaProc"), max_rows = 20)
results$print(max_rows = 20)

lazyTrace(results$draws(variables = "latent_yearly_growth[1,1]", inc_warmup = FALSE), val1 = treeData[1, dbh2 - dbh1])

lazyTrace(results$draws(variables = "latent_yearly_growth[8,2]", inc_warmup = FALSE), val1 = treeData[8, dbh3 - dbh2])

lazyTrace(results$draws(variables = "latent_yearly_growth[8,15]", inc_warmup = FALSE), val1 = treeData[8, dbh16 - dbh15])

lazyTrace(results$draws(variables = "sigmaProc", inc_warmup = FALSE), val1 = sigmaProc)

