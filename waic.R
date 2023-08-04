
#### Aim of prog: Compare the SSM approach to the classic approach using the PSIS-LOO CV approach and sensitivity analysis
## Comments
# The theory of WAIC and PSIS-LOO CV is described in Vehtari.2017:
#	Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC
#	DOI: 10.1007/s11222-016-9696-4
#	https://doi.org/10.1007/s11222-016-9696-4
#
#	Note that I stupidly reprogrammed the function waic in an older version and found exactly the same results! This version uses the package loo
#
# The theory of sensitivity analysis is described simply in Puy.2021:
#	sensobol: an R package to compute variance-based sensitivity indices
#	DOI: 10.48550/arxiv.2101.10103
#	https://arxiv.org/abs/2101.10103
#
#	Note that for a better understanding, the book Sensitivity Analysis: Primer (Saltelli, 2008) is a good option.
#	"Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index" (Saltell.2010) gives an overview 
#	of Total Sensitivity Indices.
#	DOI: 10.1016/j.cpc.2009.09.018
#	https://www.sciencedirect.com/science/article/abs/pii/S0010465509003087?via%3Dihub

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(sensobol)
library(stringi)
library(foreach)

#### Tool functions
## Source functions
source("toolFunctions.R")

## Function to compute the indices of dendro data
current_dbh = function(dendro)
{
	ls_plots = unique(dendro[, plot_id])
	count_plot = 0
	ls_trees = unique(dendro[, tree_id])

	dendro[, current_dbh := starting_dbh]
	for (current_plot in ls_plots)
	{
		for (current_tree in ls_trees)
		{
			ls_years = dendro[.(current_plot, current_tree)][2:.N, year] # Except first year
			for (current_year in ls_years)
			{
				dendro[.(current_plot, current_tree, current_year), current_dbh :=
					dendro[.(current_plot, current_tree, current_year - 1), current_dbh + dbh_increment_in_mm]]
			}
			if (!dendro[.(current_plot, current_tree)][.N, all.equal(current_dbh + dbh_increment_in_mm, dbh)])
				print(paste0("Check tree <", current_tree, "> from plot <", current_plot, ">"))
		}
		count_plot = count_plot + 1
		print(paste0(round(count_plot*100/length(ls_plots), 2), "% done"))
	}
}

## Function to check the normality of the posterior distribution of a parameter, and compute mean and sd
check_normality = function(model, params)
{
	mean_sd_dt = data.table(parameter = params, mu = numeric(length(params)), sigma = numeric(length(params)),
		shapiro_pVal = numeric(length(params)))
	setkey(mean_sd_dt, parameter)

	for (param in params)
	{
		draws = model$draws(param)
		mean_sd_dt[param, c("mu", "sigma", "shapiro_pVal") := .(mean(draws), sd(draws), shapiro.test(draws)$p.value)]
	}

	mean_sd_dt[, isNormal := shapiro_pVal > 0.05]
	return (mean_sd_dt)
}

## Function to check correlations between parameters draws
correl_draws = function(model, params, threshold = 0.75)
{
	n_iter = model$metadata()$iter_sampling
	n_chains = model$num_chains()

	mat = matrix(data = NA, nrow = n_iter*n_chains, ncol = length(params))

	for (i in seq_along(params))
		mat[, i] = as.numeric(model$draws(params[i]))

	# pairs(mat)
	N = length(params)
	
	for (i in 1:(N - 1))
		for (j in (i + 1):N)
		{
			correl = cor(mat[, i], mat[, j])
			if (abs(correl) > threshold)
				print(paste("Correlation", params[i], "--", params[j], round(correl, 2)))
		}
}

#? ----------------------------------------------------------------------------------------
#* ----------------------    PART I: Compute PSIS-LOO CV and waic    ----------------------
#? ----------------------------------------------------------------------------------------
#### Load results
## Common variables
# args = c("Betula pendula", "1")			# OK, ssm better
# args = c("Fagus sylvatica", "1")			# OK, ssm better
# args = c("Picea abies", "1")				# OK, ssm better
# args = c("Pinus pinaster", "1")			# OK, ssm better
# args = c("Pinus sylvestris", "1")			# OK, ssm better
# args = c("Quercus petraea" , "1")			# OK, ssm better
args = commandArgs(trailingOnly = TRUE) # args = c("Fagus sylvatica", "1")
if (length(args) < 2)
	stop("Supply in this order the species and the run_id", call. = FALSE)

species = as.character(args[1])
run = as.integer(args[2])

tree_path = paste0("./", species, "/")
if (!dir.exists(tree_path))
	stop(paste0("Path not found for species <", species, ">."))

## Results
info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
classic = readRDS(paste0(tree_path, info_lastRun[["file"]]))

#### Load data
## Load and subset dendrochronological data
dendro = readRDS("/home/amael/project_ssm/inventories/treeRings/treeRings-climate.rds")
dendro = dendro[speciesName_sci == species]

if (dendro[, .N] == 0)
	stop("The species is not in the dendro data set")

## Add basal area (set to zero, which translates into average basal area on the real scale)
dendro[, standBasalArea := 0] #! Considered already standardised!

## Add a key (to make sure it is sort by plot_id, tree_id, year). Here, plot_id is redundant since it is included in plot_id
setkey(dendro, plot_id, tree_id, year)

## Comput current_dbh, each row of dendro should be read this way: phi_j, gamma_j, which gives phi_{j + 1} on the next row!
current_dbh(dendro)

dendro = dendro[dbh_increment_in_mm > 0]

## Load scalings
# Load scalings
# ba_scaling_ssm = readRDS(paste0(tree_path, run, "_ba_normalisation.rds")) #! To remove?
climate_scaling_ssm = readRDS(paste0(tree_path, run, "_climate_normalisation.rds"))
ph_scaling_ssm = readRDS(paste0(tree_path, run, "_ph_normalisation.rds"))

# ba_scaling_classic = readRDS(paste0(tree_path, run, "_ba_normalisation_classic.rds")) #! To remove?
climate_scaling_classic = readRDS(paste0(tree_path, run, "_climate_normalisation_classic.rds"))
ph_scaling_classic = readRDS(paste0(tree_path, run, "_ph_normalisation_classic.rds"))

# scaling_ssm = rbindlist(list(ba_scaling_ssm, climate_scaling_ssm, ph_scaling_ssm)) #! To remove?
scaling_ssm = rbindlist(list(climate_scaling_ssm, ph_scaling_ssm))
# scaling_classic = rbindlist(list(ba_scaling_classic, climate_scaling_classic, ph_scaling_classic)) #! To remove?
scaling_classic = rbindlist(list(climate_scaling_classic, ph_scaling_classic))

scaling_dt = rbindlist(list(ssm = scaling_ssm, classic = scaling_classic), idcol = "type")
scaling_dt[, variable := stri_replace_all(str = variable, replacement = "", regex = "_avg$")] # Change the names for ease of usage
scaling_dt[variable == "pr", variable := "precip"] # Change the names for ease of usage

setkey(scaling_dt, type, variable)

## Create stan data to compute waic
# SSM
stanData_ssm = readRDS(paste0(tree_path, run, "_stanData.rds"))
stanData_ssm$n_data_rw = dendro[, .N]

stanData_ssm$precip_rw = dendro[, (pr - scaling_dt[.("ssm", "precip"), mu])/scaling_dt[.("ssm", "precip"), sd]]
stanData_ssm$tas_rw = dendro[, (tas - scaling_dt[.("ssm", "tas"), mu])/scaling_dt[.("ssm", "tas"), sd]]
stanData_ssm$ph_rw = dendro[, (ph - scaling_dt[.("ssm", "ph"), mu])/scaling_dt[.("ssm", "ph"), sd]]
stanData_ssm$standBasalArea_rw = dendro[, standBasalArea]

stanData_ssm$dbh_rw = dendro[, current_dbh/stanData_ssm$sd_dbh]
stanData_ssm$ring_width = dendro[, dbh_increment_in_mm/stanData_ssm$sd_dbh]

# Classic
stanData_classic = readRDS(paste0(tree_path, run, "_stanData_classic.rds"))
stanData_classic$n_data_rw = dendro[, .N]

stanData_classic$precip_rw = dendro[, (pr - scaling_dt[.("classic", "precip"), mu])/scaling_dt[.("classic", "precip"), sd]]
stanData_classic$tas_rw = dendro[, (tas - scaling_dt[.("classic", "tas"), mu])/scaling_dt[.("classic", "tas"), sd]]
stanData_classic$ph_rw = dendro[, (ph - scaling_dt[.("classic", "ph"), mu])/scaling_dt[.("classic", "ph"), sd]]
stanData_classic$standBasalArea_rw = dendro[, standBasalArea]

stanData_classic$dbh_rw = dendro[, current_dbh/stanData_classic$sd_dbh]
stanData_classic$ring_width = dendro[, dbh_increment_in_mm/stanData_classic$sd_dbh]

#### Compute likelihood
## Common variables
model = cmdstanr::cmdstan_model("./waic.stan")
n_chains = ssm$num_chains()
n_iter = ssm$metadata()$iter_sampling

loglik_ssm = model$generate_quantities(ssm$draws(), data = stanData_ssm, parallel_chains = n_chains)
loglik_classic = model$generate_quantities(classic$draws(), data = stanData_classic, parallel_chains = n_chains)

r_eff_ssm = loo::relative_eff(loglik_ssm$draws("log_lik"), cores = 8)
loo_ssm = loo::loo(x = loglik_ssm$draws("log_lik"), r_eff = r_eff_ssm, cores = 8)
waic_ssm = loo::waic(x = loglik_ssm$draws("log_lik"), cores = 8)

r_eff_classic = loo::relative_eff(loglik_classic$draws("log_lik"), cores = 8)
loo_classic = loo::loo(x = loglik_classic$draws("log_lik"), r_eff = r_eff_classic, cores = 8)
waic_classic = loo::waic(x = loglik_classic$draws("log_lik"), cores = 8)

loo::loo_compare(list(loo_ssm, loo_classic))
loo::loo_compare(list(waic_ssm, waic_classic))

#? -----------------------------------------------------------------------------------------
#* ----------------------    PART II: Compute sensitivity analysis    ----------------------
#? -----------------------------------------------------------------------------------------
sensitivityAnalysis = function(model, dbh0, sd_dbh, env0, matrices = c("A", "B", "AB"), order = "second", N = 2^14, boot = FALSE,
	bootstrap_iter = NULL, type = "QRN")
{
	## Check dbh0 and env0
	if (any(abs(c(dbh0, env0)) > 3))
		warning("Are you sure that dbh0 and env0 are scaled?")

	## Check boot
	if (boot && is.null(bootstrap_iter))
	{
		warning("Boot is true, but no boostrap_iter provided! A value of 3000 was assigned by default")
		bootstrap_iter = 3000
	}

	if (!boot && !is.null(bootstrap_iter))
	{
		warning("Boot is false, so the provided value for bootstrap_iter is disregarded")
		bootstrap_iter = NULL
	}
	
	## Prepare the sampling matrices
	# Common variables
	params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2",
		"competition_slope")

	explanatory_vars = "dbh"

	# Check normality
	mean_sd_dt = check_normality(model = model, params = params)

	if (any(mean_sd_dt[, !isNormal]))
	{
		print("The following are detected to be non-normal:")
		print(mean_sd_dt[(!isNormal), parameter])
	}

	## Create matrices
	sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = c(params, explanatory_vars))

	# Rescale the Sobol matrix
	# --- Parameters
	for (param in params)
		sobol_mat[, param] = qnorm(sobol_mat[, param], mean_sd_dt[param, mu], mean_sd_dt[param, sigma])

	# --- Diameters
	sobol_mat[, "dbh"] = qnorm(p = sobol_mat[, "dbh"], mean = dbh0, sd = 3/sd_dbh)

	if (any(is.na(sobol_mat)))
		stop(paste("Number of NAs in Sobol's matrix:", sum(is.na(sobol_mat))))

	## Run sensitivity analysis
	# Model output
	y = numeric(nrow(sobol_mat))
	for (i in seq_len(nrow(sobol_mat)))
	{
		params_vec = c(sobol_mat[i, "averageGrowth"],
			sobol_mat[i, "dbh_slope"],
			sobol_mat[i, "dbh_slope2"],
			sobol_mat[i, "pr_slope"],
			sobol_mat[i, "tas_slope"],
			sobol_mat[i, "ph_slope"],
			sobol_mat[i, "competition_slope"],
			sobol_mat[i, "pr_slope2"],
			sobol_mat[i, "tas_slope2"],
			sobol_mat[i, "ph_slope2"])
		
		dbh = sobol_mat[i, "dbh"]
		
		y[i] = growth_fct_meanlog(dbh = dbh, pr = env0["precip"], tas = env0["tas"], ph = env0["ph"], basalArea = env0["basalArea"],
			params = params_vec, sd_dbh = sd_dbh, standardised_variables = TRUE)

		if (i %% 20000 == 0)
			print(paste0(round(100*i/nrow(sobol_mat), 2), "% done"))
	}
	print("100% done")

	# Sobol indices
	ind = sobol_indices(matrices = matrices, Y = y, N = N, params = c(params, "dbh"), first = "saltelli", total = "jansen",
		boot = boot, order = order, R = bootstrap_iter, parallel = "multicore", type = "percent", conf = 0.95, ncpus = 8)

	return(ind)
}

sa_ssm = sensitivityAnalysis(model = ssm, dbh0 = 1.2, sd_dbh = stanData_ssm$sd_dbh, order = "first",
	env0 = c(precip = 0.8, tas = -0.5, ph = 0.21, basalArea = 0.58), N = 2^16)$results
sa_ssm = data.table::dcast(data = sa_ssm, formula = parameters ~ sensitivity, value.var = "original")
sa_ssm[, sum(Si)]
sa_ssm = sa_ssm[order(-Si), .SD]

# total_ssm = sa_ssm$results[sensitivity == "Ti"]
# total_ssm = total_ssm[order(-original), .SD, by = sensitivity]
# saveRDS(sa_ssm, paste0(tree_path, run, "_sa_ssm"))

sa_classic = sensitivityAnalysis(model = classic, dbh0 = 1.2, sd_dbh = stanData_classic$sd_dbh, order = "first",
	env0 = c(precip = 0.8, tas = -0.5, ph = 0.21, basalArea = 0.58), N = 2^16)$results
sa_classic = data.table::dcast(data = sa_classic, formula = parameters ~ sensitivity, value.var = "original")
sa_classic[, sum(Si)]
sa_classic = sa_classic[order(-Si), .SD]

# total_classic = sa_classic$results[sensitivity == "Ti"]
# total_classic = total_classic[order(-original), .SD, by = sensitivity]
# saveRDS(sa_classic, paste0(tree_path, run, "_sa_classic"))

comp = merge.data.table(sa_ssm, sa_classic, by = "parameters", suffixes = c("_ssm", "_classic"))
comp[, diff_Si := 100*(Si_ssm - Si_classic)/Si_classic]
comp[, diff_Ti := 100*(Ti_ssm - Ti_classic)/Ti_classic]

comp[, .(parameters, diff_Si, diff_Ti)]


#? ------------------------------------------------------------------------------------------------------------------------
#* ----------------------    PART III: Compute sensitivity analysis with respect to climate input    ----------------------
#? ------------------------------------------------------------------------------------------------------------------------

fraction_range_climate = function(stanData, model_type, shouldScale = FALSE, frac = 0.01, scaling_dt = NULL)
{
	if (shouldScale && is.null(scaling_dt))
		stop("Scaling demanded, but no scaling provided!")
	
	range_dt = data.table(variable = c("precip", "tas", "ph"), q05 = numeric(3), q95 = numeric(3))
	setkey(range_dt, variable)
	for (currentVar in range_dt[, variable])
	{
		clim = stanData_ssm[[currentVar]]
		if (shouldScale)
			clim = (clim - scaling_dt[.(model_type, currentVar), mu])/scaling_dt[.(model_type, currentVar), sd]
		range_dt[currentVar, c("q05", "q95") := as.list(quantile(clim, c(0.05, 0.95)))]
	}

	range_dt[, range := q95 - q05]
	range_dt[, sd_sa := frac*range]

	return(range)
}

