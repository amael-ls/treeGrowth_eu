
#### Aim of prog: Replicate the small simulation done in Aubry-Kientz et. al (2017) to verify if Eitzel 2013 can be applied for Δt >= 5

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

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

## Data
treeData = data.table(tree_id = 1:n, dbh1 = init_dbh)
treeData[, paste0("dbh", 2:n_years) := 0]

for (i in 2:n_years)
{
	current_dbh = paste0("dbh", i)
	previous_dbh = paste0("dbh", i - 1)

	treeData[, c(current_dbh) := (1 + beta1)*treeData[, ..previous_dbh] + beta0 + beta2*temperature[i]]
}

#### Stan model
## Compile model
model = cmdstan_model("Aubry-Kientz2017.stan")

## Stan data
# Common variables
freq_obs = 2

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

## Run model
n_chains = 1
model = cmdstan_model("Aubry-Kientz2017.stan")
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 1500, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 10, adapt_delta = 0.85)

results$cmdstan_diagnose()

results$print(c("lp__", "beta0", "beta1", "beta2", "sigmaProc"), max_rows = 20)
