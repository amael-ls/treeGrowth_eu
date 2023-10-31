
#### Aim of prog: Analysing results (check-up residuals, plots)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(reticulate)
library(posterior)
library(cmdstanr)
library(stringi)
library(png)

#### Tool functions
source("./toolFunctions.R")

#### Common variables
## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

n_runs = 4 # Number of runs used in growth_subsample.R
threshold_indiv = 12000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2023/01/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_de-fr-sw_12000_main.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

error_ls = vector(mode = "list", length = infoSpecies[, .N])
names(error_ls) = infoSpecies[, speciesName_sci]

correl_ls = vector(mode = "list", length = infoSpecies[, .N])
names(correl_ls) = infoSpecies[, speciesName_sci]

posterior_ls = vector(mode = "list", length = infoSpecies[, .N])
names(posterior_ls) = infoSpecies[, speciesName_sci]

ls_files = character(length = infoSpecies[, .N])
names(ls_files) = infoSpecies[, speciesName_sci]

params_dt = data.table(parameters = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2",
	"tas_slope", "tas_slope2", "ph_slope", "ph_slope2", "competition_slope"),
	priors = c(dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm),
	arg1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
	arg2 = c(20, 20, 20, 20, 20, 20, 20, 20, 20, 20),
	title = c("Average growth (mean)", "Dbh slope", "Dbh slope (quadratic term)", "Precipitation slope",
		"Precipitation slope (quadratic term)", "Temperature slope", "Temperature slope (quadratic term)", "Soil acidity slope (pH)",
		"Soil acidity slope (pH, quadratic term)", "Competition slope"),
	expand_bounds = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

setkey(params_dt, parameters)

#### For loop on processed species, to plot posteriors of the main parameters (errors, intercept, and slopes)
for (species in infoSpecies[, speciesName_sci])
{
	multi = infoSpecies[species, multiRun]
	ls_nfi = unlist(stri_split(infoSpecies[species, ls_nfi], regex = ", "))
	summary_dt = centralised_fct(species, multi, n_runs, ls_nfi, params_dt, run = if (multi) NULL else 1, simulatePosterior = FALSE)
	error_ls[[species]] = summary_dt[["error_dt"]]
	correl_ls[[species]] = summary_dt[["correl_energy"]]
	posterior_ls[[species]] = summary_dt[["posteriorSim"]]
	ls_files[[species]] = summary_dt[["fileResults"]]
}

error_dt = rbindlist(error_ls, idcol = "speciesName_sci")
correl_dt = rbindlist(correl_ls, idcol = "speciesName_sci")

saveRDS(error_dt, "./error_species.rds")
saveRDS(correl_dt, "./correlation_energy_species.rds")
saveRDS(posterior_ls, "./posterior_ls.rds")

plot_correl_error(error_dt, correl_dt, threshold_correl = 0.2, rm_correl = "lp__")
