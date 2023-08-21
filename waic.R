
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

## Function to compute the climate range from stan data
fraction_range_climate = function(stanData, model_type, shouldScale = FALSE, frac = 0.01, scaling_dt = NULL)
{
	if (shouldScale && is.null(scaling_dt))
		stop("Scaling demanded, but no scaling provided!")
	
	range_dt = data.table(variable = c("precip", "tas", "ph", "standBasalArea"), q05 = numeric(4), q95 = numeric(4), mu = numeric(4))
	setkey(range_dt, variable)
	for (currentVar in range_dt[, variable])
	{
		clim = stanData[[currentVar]]
		if (shouldScale)
			clim = (clim - scaling_dt[.(model_type, currentVar), mu])/scaling_dt[.(model_type, currentVar), sd]
		range_dt[currentVar, c("q05", "q95", "mu") := as.list(c(quantile(clim, c(0.05, 0.95)), mean(clim)))]
	}

	range_dt[, range := q95 - q05]
	range_dt[, sd_sa := frac*range]

	return(range_dt)
}

## Function to compute sensitivity of growth with respect to uncertainty in the parameters and data
sensitivityAnalysis = function(model, dbh0, sd_dbh, env0, clim_dt, matrices = c("A", "B", "AB"), order = "first", N = 2^14, boot = FALSE,
	bootstrap_iter = NULL, type = "QRN", first = "saltelli", total = "jansen")
{
	## Check dbh0 and env0
	if (abs(dbh0) > 5)
		warning("Are you sure that dbh0 is scaled?")

	if ((env0["precip"] < clim_dt["precip", q05]) || (env0["tas"] < clim_dt["tas", q05]) || (env0["ph"] < clim_dt["ph", q05]))
		warning("Are you sure that env0 is scaled? There are low values below the 5 quantile")

	if ((env0["precip"] > clim_dt["precip", q95]) || (env0["tas"] > clim_dt["tas", q95]) || (env0["ph"] > clim_dt["ph", q95]))
		warning("Are you sure that env0 is scaled? There are large values beyond the 95 quantile")

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
	
	explanatory_vars = c("dbh", "precip", "tas", "ph")

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

	# --- Explanatory variables
	for (currentVar in explanatory_vars)
	{
		if (currentVar == "dbh")
		{
			sobol_mat[, "dbh"] = qnorm(p = sobol_mat[, "dbh"], mean = dbh0, sd = 3/sd_dbh)
		} else {
			sobol_mat[, currentVar] = qnorm(p = sobol_mat[, currentVar], mean = env0[currentVar], sd = clim_dt[currentVar, sd_sa])
		}
	}
	
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
		
		y[i] = growth_fct_meanlog(dbh = sobol_mat[i, "dbh"], pr = sobol_mat[i, "precip"], tas = sobol_mat[i, "tas"],
			ph = sobol_mat[i, "ph"], basalArea = env0["basalArea"],
			params = params_vec, sd_dbh = sd_dbh, standardised_variables = TRUE)

		if (i %% 20000 == 0)
			print(paste0(round(100*i/nrow(sobol_mat), 2), "% done"))
	}
	print("100% done")

	# Sobol indices
	ind = sobol_indices(matrices = matrices, Y = y, N = N, params = c(params, explanatory_vars), first = first, total = total,
		boot = boot, order = order, R = bootstrap_iter, parallel = "multicore", type = "percent", conf = 0.95, ncpus = 8)

	return(list(sa = ind, var_y = var(y)))
}

## Function to compute sensitivity of growth with respect to uncertainty in the data only
sensitivityAnalysis_data = function(model, dbh0, sd_dbh, clim_dt, n_param, env0 = NULL, matrices = c("A", "B", "AB"), order = "first", N = 2^14,
	type = "QRN", first = "saltelli", total = "jansen", seed = NULL)
{
	## Check dbh0 and env0
	if (any(abs(dbh0) > 5))
		warning("Are you sure that dbh0 is scaled?")

	if (!is.null(env0))
	{
		if ((env0["precip"] < clim_dt["precip", q05]) || (env0["tas"] < clim_dt["tas", q05]) || (env0["ph"] < clim_dt["ph", q05]))
			warning("Are you sure that env0 is scaled? There are low values below the 5 quantile")

		if ((env0["precip"] > clim_dt["precip", q95]) || (env0["tas"] > clim_dt["tas", q95]) || (env0["ph"] > clim_dt["ph", q95]))
			warning("Are you sure that env0 is scaled? There are large values beyond the 95 quantile")
	}

	if ((length(dbh0) != 1) && (length(dbh0) != 2))
	{
		warning("Only the first two values of dbh0 are used")
		dbh0 = dbh0[1:2]
	}

	if (length(dbh0) == 2)
	{
		if (dbh0[1] > dbh0[2])
		{
			dbh0[1] = dbh0[2] - dbh0[1] + (dbh0[2] = dbh0[1])
			warning("The values for dbh0 were swaped")
		}
	}
	
	## Prepare the sampling matrices
	# Common variables
	params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2",
		"competition_slope")

	if (!is.null(seed))
		set.seed(seed)
	
	params_values = getParams(model_cmdstan = model, params_names = params, type = "all")
	
	n_iter = model$metadata()$iter_sampling
	n_chains = model$num_chains()
	n_tot = n_iter*n_chains
	
	sample_ind = sample(x = 1:n_tot, size = n_param, replace = FALSE)

	explanatory_vars = c("dbh", "precip", "tas", "ph", "standBasalArea")

	## Create matrices
	sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = c(explanatory_vars))

	# Rescale the Sobol matrix
	# --- Explanatory variables
	for (currentVar in explanatory_vars)
	{
		if (currentVar == "dbh")
		{
			if (length(dbh0) == 1)
			{
				sobol_mat[, "dbh"] = qnorm(p = sobol_mat[, "dbh"], mean = dbh0, sd = 3/sd_dbh)
			} else {
				sobol_mat[, "dbh"] = qunif(p = sobol_mat[, "dbh"], min = dbh0[1], max = dbh0[2])
			}
		} else {
			if (!is.null(env0))
			{
				sobol_mat[, currentVar] = qnorm(p = sobol_mat[, currentVar], mean = env0[currentVar], sd = clim_dt[currentVar, sd_sa])
			} else {
				sobol_mat[, currentVar] = qunif(p = sobol_mat[, currentVar], min = clim_dt[currentVar, q05], max = clim_dt[currentVar, q95])
			}
		}
	}
	
	if (any(is.na(sobol_mat)))
		stop(paste("Number of NAs in Sobol's matrix:", sum(is.na(sobol_mat))))

	## Run sensitivity analysis
	ind = vector(mode = "list", length = n_param)
	var_vec_y = numeric(length = n_param)
	range_dt = data.table(run = seq_along(sample_ind), min_y = -Inf, max_y = +Inf)
	
	freq_print = 1
	if (n_param > 19)
		freq_print = round(5*n_param/100)

	for (j in seq_along(sample_ind))
	{
		current_ind = sample_ind[j]
		params_vec = c(averageGrowth = params_values[, , "averageGrowth"][current_ind],
			dbh_slope = params_values[, , "dbh_slope"][current_ind],
			dbh_slope2 = params_values[, , "dbh_slope2"][current_ind],
			pr_slope = params_values[, , "pr_slope"][current_ind],
			tas_slope = params_values[, , "tas_slope"][current_ind],
			ph_slope = params_values[, , "ph_slope"][current_ind],
			competition_slope = params_values[, , "competition_slope"][current_ind],
			pr_slope2 = params_values[, , "pr_slope2"][current_ind],
			tas_slope2 = params_values[, , "tas_slope2"][current_ind],
			ph_slope2 = params_values[, , "ph_slope2"][current_ind])
		
		# --- Model output
		y = numeric(nrow(sobol_mat))
		for (i in seq_len(nrow(sobol_mat)))
		{	
			y[i] = growth_fct_meanlog(dbh = sobol_mat[i, "dbh"], pr = sobol_mat[i, "precip"], tas = sobol_mat[i, "tas"],
				ph = sobol_mat[i, "ph"], basalArea = sobol_mat[i, "standBasalArea"],
				params = params_vec, sd_dbh = sd_dbh, standardised_variables = TRUE)

		}
		var_vec_y[j] = var(y)
		range_dt[j, c("min_y", "max_y") := .(min(y), max(y))]

		# --- Sobol indices
		ind[[j]] = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total,
			order = order)$results

		if (j %% freq_print == 0)
			print(paste0(round(100*j/n_param), "% done"))
	}
	print("100% done")

	ind = rbindlist(l = ind, idcol = "run")
	if (any(ind[, original] < 0))
			warning(paste("There are negatives Si, with min value:", round(ind[, min(original)], 7)))

	return(list(sa = ind, var_y = var_vec_y, n_param = n_param, range_dt = range_dt))
}

#? ----------------------------------------------------------------------------------------
#* ----------------------    PART I: Compute PSIS-LOO CV and waic    ----------------------
#? ----------------------------------------------------------------------------------------
#### Load results
## Common variables
# args = c("Betula pendula", "1", "dbhm", "precipm", "tasm")				# OK, ssm better
# args = c("Betula pendula", "1")											# OK, ssm better
# args = c("Fagus sylvatica", "1", "dbhm", "precipm", "tasm")				# OK, ssm better
# args = c("Fagus sylvatica", "1")											# OK, ssm better
# args = c("Picea abies", "1", "dbhm", "precipm", "tasm")					# OK, ssm better
# args = c("Picea abies", "1")												# OK, ssm better
# args = c("Pinus pinaster", "1", "dbhm", "precipm", "tasm")				# OK, ssm better
# args = c("Pinus pinaster", "1")											# OK, ssm better
# args = c("Pinus sylvestris", "1", "dbhm", "precipm", "tasm")				# OK, ssm better
# args = c("Pinus sylvestris", "1")											# OK, ssm better
# args = c("Quercus petraea" , "1", "dbhm", "precipm", "tasm")				# OK, ssm better
# args = c("Quercus petraea" , "1")											# OK, ssm better
args = commandArgs(trailingOnly = TRUE)
for (i in seq_along(args))
	print(paste0("Arg ", i, ": <", args[i], ">"))

args[1] = paste(args[1], args[2])
args = args[-2]

# if (length(args) != 5)
# 	stop("Supply in this order the species, the run_id, and the quantiles for dbh, tas, and precip", call. = FALSE)
if (length(args) != 2)
	stop("Supply in this order the species, the run_id, and the quantiles for dbh, tas, and precip", call. = FALSE)


(species = as.character(args[1]))
(run = as.integer(args[2]))
# (dbh_quantile_opt = as.character(args[3]))
# (pr_quantile_opt = as.character(args[4]))
# (tas_quantile_opt = as.character(args[5]))

# ls_opt_Qt = paste0(rep(c("dbh", "tas", "precip"), each = 2), c("m", "M"))

# if (!all(c(dbh_quantile_opt, pr_quantile_opt, tas_quantile_opt) %in% ls_opt_Qt))
# 	stop(paste("Only the following options can be used:", paste(ls_opt_Qt, collapse = ", ")))

# opt_code = paste0(stri_sub(str = c(dbh_quantile_opt, pr_quantile_opt, tas_quantile_opt), from = -1), collapse = "_")

# qt_opt = c("dbh" = ifelse(dbh_quantile_opt == "dbhm", "5%", "95%"),
# 	"precip" = ifelse(pr_quantile_opt == "precipm", "q05", "q95"),
# 	"tas" = ifelse(tas_quantile_opt == "tasm", "q05", "q95"))

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
# Load climate and ph scalings (dbh and basal area useless here, as BA already standardised, and the good sd_dbh is in stanData)
climate_scaling_ssm = readRDS(paste0(tree_path, run, "_climate_normalisation.rds"))
ph_scaling_ssm = readRDS(paste0(tree_path, run, "_ph_normalisation.rds"))
ba_scaling_ssm = readRDS(paste0(tree_path, run, "_ba_normalisation.rds"))

climate_scaling_classic = readRDS(paste0(tree_path, run, "_climate_normalisation_classic.rds"))
ph_scaling_classic = readRDS(paste0(tree_path, run, "_ph_normalisation_classic.rds"))
ba_scaling_classic = readRDS(paste0(tree_path, run, "_ba_normalisation_classic.rds"))

scaling_ssm = rbindlist(list(climate_scaling_ssm, ph_scaling_ssm, ba_scaling_ssm))
scaling_classic = rbindlist(list(climate_scaling_classic, ph_scaling_classic, ba_scaling_classic))

scaling_dt = rbindlist(list(ssm = scaling_ssm, classic = scaling_classic), idcol = "type")
scaling_dt[, variable := stri_replace_all(str = variable, replacement = "", regex = "_avg$")] # Change the names for ease of usage
scaling_dt[variable == "pr", variable := "precip"] # Change the names for ease of usage
scaling_dt[variable == "standBasalArea_interp", variable := "standBasalArea"] # Change the names for ease of usage

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
#### Climate range
clim_dt_ssm = fraction_range_climate(stanData = stanData_ssm, model_type = "ssm", shouldScale = TRUE, scaling_dt = scaling_dt)
clim_dt_classic = fraction_range_climate(stanData = stanData_ssm, model_type = "classic", shouldScale = TRUE, scaling_dt = scaling_dt)
#! I want to have the same climate range, so even for classic, I use ssm! Actually this line is useless, it is just for backward compatibility

# clim_dt_ssm = merge.data.table(x = clim_dt_ssm, y = scaling_dt[("ssm")], by = "variable")
# clim_dt_classic = merge.data.table(x = clim_dt_classic, y = scaling_dt[("classic")], by = "variable")

# clim_dt_ssm[, mu.x := mu.x*sd + mu.y]
# clim_dt_ssm[, q05 := q05*sd + mu.y]
# clim_dt_ssm[, q95 := q95*sd + mu.y]

# clim_dt_classic[, mu.x := mu.x*sd + mu.y]
# clim_dt_classic[, q05 := q05*sd + mu.y]
# clim_dt_classic[, q95 := q95*sd + mu.y]

# dbh0_ssm = unname(quantile(stanData_ssm$dbh_init, c(0.05, 0.95))[qt_opt["dbh"]])/stanData_ssm$sd_dbh
# dbh0_classic = unname(quantile(stanData_classic$dbh_init, c(0.05, 0.95))[qt_opt["dbh"]])/stanData_classic$sd_dbh

# #### Sensitivity of growth with respect to parameters and data
# ## SSM
# sa_ssm = sensitivityAnalysis(model = ssm, dbh0 = dbh0_ssm, sd_dbh = stanData_ssm$sd_dbh,
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0),
# 	clim_dt = clim_dt_ssm, N = 2^19)
# saveRDS(sa_ssm, paste0(tree_path, "sa_ssm_", opt_code, ".rds"))

# sa_ssm = data.table::dcast(data = sa_ssm$results, formula = parameters ~ sensitivity, value.var = "original")
# sa_ssm[, sum(Si)]
# sa_ssm = sa_ssm[order(-Si), .SD]

# if (any(sa_ssm[, Si] < 0))
# 	warning(paste("There are negatives Si for SSM, with min value:", round(sa_ssm[, min(Si)], 7)))

# ## Classic
# sa_classic = sensitivityAnalysis(model = classic, dbh0 = dbh0_classic, sd_dbh = stanData_classic$sd_dbh,
# 	env0 = c(precip = clim_dt_classic["precip", get(qt_opt["precip"])], tas = clim_dt_classic["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0),
# 	clim_dt = clim_dt_classic, N = 2^19)
# saveRDS(sa_classic, paste0(tree_path, "sa_classic_", opt_code, ".rds"))

# sa_classic = data.table::dcast(data = sa_classic$results, formula = parameters ~ sensitivity, value.var = "original")
# sa_classic[, sum(Si)]
# sa_classic = sa_classic[order(-Si), .SD]

# if (any(sa_classic[, Si] < 0))
# 	warning(paste("There are negatives Si for classic, with min value:", round(sa_classic[, min(Si)], 7)))

# [1] "The following are detected to be non-normal:"
# [1] "averageGrowth" "dbh_slope"     "dbh_slope2"   
# [4] "ph_slope"      "ph_slope2"     "pr_slope2"    
# [7] "tas_slope2" 

#### Sensitivity of growth with respect to data only
## SSM
sa_ssm_data = sensitivityAnalysis_data(model = ssm, dbh0 = quantile(stanData_ssm$dbh_init, c(0.05, 0.95))/stanData_ssm$sd_dbh,
	sd_dbh = stanData_ssm$sd_dbh, n_param = 500, N = 2^16, clim_dt = clim_dt_ssm, seed = 123)
# sa_ssm_data = sensitivityAnalysis_data(model = ssm, dbh0 = dbh0_ssm, sd_dbh = stanData_ssm$sd_dbh, n_param = 1e3, N = 2^16,
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0),
# 	clim_dt = clim_dt_ssm, seed = 123)

saveRDS(sa_ssm_data, paste0(tree_path, "sa_ssm_data.rds"))
# saveRDS(sa_ssm_data, paste0(tree_path, "sa_ssm_", opt_code, "_data.rds"))

## Classic
sa_classic_data = sensitivityAnalysis_data(model = classic, dbh0 = quantile(stanData_ssm$dbh_init, c(0.05, 0.95))/stanData_ssm$sd_dbh,
	sd_dbh = stanData_classic$sd_dbh, n_param = 500, N = 2^16, clim_dt = clim_dt_classic, seed = 123)
# sa_classic_data = sensitivityAnalysis_data(model = classic, dbh0 = dbh0_classic, sd_dbh = stanData_classic$sd_dbh, n_param = 1e3, N = 2^16,
# 	env0 = c(precip = clim_dt_classic["precip", get(qt_opt["precip"])], tas = clim_dt_classic["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0),
# 	clim_dt = clim_dt_classic, seed = 123)

saveRDS(sa_classic_data, paste0(tree_path, "sa_classic_data.rds"))
# saveRDS(sa_classic_data, paste0(tree_path, "sa_classic_", opt_code, "_data.rds"))

# setkey(sa_ssm_data$sa, parameters, sensitivity)
# setkey(sa_classic_data$sa, parameters, sensitivity)

# plot_sa(sa_ssm_data$sa, sa_classic_data$sa, dataNames = c("dbh", "precip", "tas", "ph"), plot_opt = "separate", type = "Si", n = 2048)


print(sa_ssm_data$sa[sensitivity == "Si", sum(original), by = run][, range(V1)])
print(sa_classic_data$sa[sensitivity == "Si", sum(original), by = run][, range(V1)])



# ## Comparison
# comp = merge.data.table(sa_ssm, sa_classic, by = "parameters", suffixes = c("_ssm", "_classic"))
# comp[, diff_Si := 100*(Si_ssm - Si_classic)/Si_classic]
# comp[, diff_Ti := 100*(Ti_ssm - Ti_classic)/Ti_classic]

# comp[, .(parameters, diff_Si, diff_Ti)]
# comp = comp[order(-Si_ssm), .SD]
# comp[, .(parameters, diff_Si)]



# ## Check second order
# sa_ssm_2 = sensitivityAnalysis(model = ssm, dbh0 = dbh0_ssm, sd_dbh = stanData_ssm$sd_dbh, order = "second",
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0), clim_dt = clim_dt_ssm, N = 2^19)$results
# sa_ssm_2[, sum(original), by = sensitivity] #! Shows that there is no interaction, which was expected! 
# sa_ssm_2[parameters == "dbh_slope.dbh"]
# sa_ssm_2[parameters == "dbh"]
# sa_ssm_2[sensitivity == "Si", sum(original)] + sa_ssm_2[sensitivity == "Sij", sum(original)]

# ####! CRASH TEST ZONE
# params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
# 	"ph_slope", "ph_slope2", "competition_slope")
# explanatory_vars = c("dbh", "precip", "tas", "ph")

# ## Function to create the "product names", such as dbh with dbh_slope, tas with tas_slope and tas_slope2, etc...
# product_names = function(search_in, to_search)
# {
# 	results = c()
# 	for (reg in to_search)
# 	{
# 		save_reg = reg
# 		if (reg == "precip")
# 			reg = "pr"
# 		results = append(x = results, values = paste0(search_in[stri_detect(str = search_in, regex = reg)], ".", save_reg))
# 	}

# 	return (results)
# }

# cross_names = product_names(search_in = params, to_search = explanatory_vars)

# ## SSM, 2nd order
# sa_ssm_2 = sensitivityAnalysis(model = ssm, dbh0 = dbh0_ssm, sd_dbh = stanData_ssm$sd_dbh, order = "second",
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0), clim_dt = clim_dt_ssm, N = 2^19)

# saveRDS(sa_ssm_2, paste0(tree_path, "sa_ssm_2.rds"))

# sa_ssm_2 = sa_ssm_2$results
# if (any(sa_ssm_2[, original] < 0))
# 	warning(paste("There are negatives values for ssm second order, with a minimum of", round(sa_ssm_2[, min(original)], 5)))

# sa_ssm_2[, sum(original), by = sensitivity]
# sa_ssm_2[sensitivity == "Si", sum(original)] + sa_ssm_2[sensitivity == "Sij", sum(original)]

# sa_ssm_2[parameters %in% cross_names]
# sa_ssm_2[sensitivity == "Sij"][!(parameters %in% cross_names), sum(original)]

# Si_ssm = data.table::dcast(data = sa_ssm_2[sensitivity %in% c("Si", "Ti")], formula = parameters ~ sensitivity, value.var = "original")
# Sij_ssm = data.table::dcast(data = sa_ssm_2[sensitivity == "Sij"], formula = parameters ~ sensitivity, value.var = "original")
# Sij_ssm_cross = data.table::dcast(data = sa_ssm_2[parameters %in% cross_names], formula = parameters ~ sensitivity, value.var = "original")

# Si_ssm = Si_ssm[order(-Si), .SD]
# Sij_ssm = Sij_ssm[order(-Sij), .SD]
# Sij_ssm_cross = Sij_ssm_cross[order(-Sij), .SD]


# comp_1_2_order = merge.data.table(x = Si_ssm, y = sa_ssm, by = "parameters", suffixes = c("_2", "_1"))
# comp_1_2_order[, Si_ratio_percent := 100*Si_2/Si_1]


# ## By themes (cf p. 37, Saltelli 2008)
# ssm_from_params = sa_ssm[parameters %in% params, sum(Si)]
# ssm_from_data = sa_ssm[parameters %in% explanatory_vars, sum(Si)]

# classic_from_params = sa_classic[parameters %in% params, sum(Si)]
# classic_from_data = sa_classic[parameters %in% explanatory_vars, sum(Si)]

# 100*(ssm_from_params - classic_from_params)/classic_from_params
# 100*(ssm_from_data - classic_from_data)/classic_from_data

# ####! END CRASH TEST ZONE


# ####! CRASH TEST ZONE II: test posterior SA(G, data)
# ## Function to compute sensitivity of growth with respect to uncertainty in the data only
# sensitivityAnalysis_data = function(model, dbh0, sd_dbh, env0, clim_dt, n_param, matrices = c("A", "B", "AB"), order = "first", N = 2^14,
# 	type = "QRN", first = "saltelli", total = "jansen", seed = NULL)
# {
# 	## Check dbh0 and env0
# 	if (abs(dbh0) > 5)
# 		warning("Are you sure that dbh0 is scaled?")

# 	if ((env0["precip"] < clim_dt["precip", q05]) || (env0["tas"] < clim_dt["tas", q05]) || (env0["ph"] < clim_dt["ph", q05]))
# 		warning("Are you sure that env0 is scaled? There are low values below the 5 quantile")

# 	if ((env0["precip"] > clim_dt["precip", q95]) || (env0["tas"] > clim_dt["tas", q95]) || (env0["ph"] > clim_dt["ph", q95]))
# 		warning("Are you sure that env0 is scaled? There are large values beyond the 95 quantile")
	
# 	## Prepare the sampling matrices
# 	# Common variables
# 	params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2",
# 		"competition_slope")

# 	if (!is.null(seed))
# 		set.seed(seed)
	
# 	params_values = getParams(model_cmdstan = model, params_names = params, type = "all")
	
# 	n_iter = model$metadata()$iter_sampling
# 	n_chains = model$num_chains()
# 	n_tot = n_iter*n_chains
	
# 	sample_ind = sample(x = 1:n_tot, size = n_param, replace = FALSE)

# 	explanatory_vars = c("dbh", "precip", "tas", "ph")

# 	## Create matrices
# 	sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = c(explanatory_vars))

# 	# Rescale the Sobol matrix
# 	# --- Explanatory variables
# 	for (currentVar in explanatory_vars)
# 	{
# 		if (currentVar == "dbh")
# 		{
# 			sobol_mat[, "dbh"] = qnorm(p = sobol_mat[, "dbh"], mean = dbh0, sd = 3/sd_dbh)
# 		} else {
# 			sobol_mat[, currentVar] = qnorm(p = sobol_mat[, currentVar], mean = env0[currentVar], sd = clim_dt[currentVar, sd_sa])
# 		}
# 	}
	
# 	if (any(is.na(sobol_mat)))
# 		stop(paste("Number of NAs in Sobol's matrix:", sum(is.na(sobol_mat))))

# 	## Run sensitivity analysis
# 	ind = vector(mode = "list", length = n_param)
# 	for (j in seq_along(sample_ind))
# 	{
# 		current_ind = sample_ind[j]
# 		params_vec = c(averageGrowth = params_values[, , "averageGrowth"][current_ind],
# 			dbh_slope = params_values[, , "dbh_slope"][current_ind],
# 			dbh_slope2 = params_values[, , "dbh_slope2"][current_ind],
# 			pr_slope = params_values[, , "pr_slope"][current_ind],
# 			tas_slope = params_values[, , "tas_slope"][current_ind],
# 			ph_slope = params_values[, , "ph_slope"][current_ind],
# 			competition_slope = params_values[, , "competition_slope"][current_ind],
# 			pr_slope2 = params_values[, , "pr_slope2"][current_ind],
# 			tas_slope2 = params_values[, , "tas_slope2"][current_ind],
# 			ph_slope2 = params_values[, , "ph_slope2"][current_ind])
		
# 		# --- Model output
# 		y = numeric(nrow(sobol_mat))
# 		for (i in seq_len(nrow(sobol_mat)))
# 		{	
# 			y[i] = growth_fct_meanlog(dbh = sobol_mat[i, "dbh"], pr = sobol_mat[i, "precip"], tas = sobol_mat[i, "tas"],
# 				ph = sobol_mat[i, "ph"], basalArea = env0["basalArea"],
# 				params = params_vec, sd_dbh = sd_dbh, standardised_variables = TRUE)

# 		}

# 		# --- Sobol indices
# 		ind[[j]] = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total,
# 			order = order)$results

# 		if (j %% 5 == 0)
# 			print(paste0(round(100*j/n_param, 2), "% done"))
# 	}
# 	print("100% done")

# 	ind = rbindlist(l = ind, idcol = "run")
# 	if (any(ind[, original] < 0))
# 			warning("There are negatives Si")

# 	return(ind)
# }

# sa_ssm_data = sensitivityAnalysis_data(model = ssm, dbh0 = dbh0_ssm, sd_dbh = stanData_ssm$sd_dbh, n_param = 1e3, N = 2^14,
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0), clim_dt = clim_dt_ssm, seed = 123)

# saveRDS(sa_ssm_data, paste0(tree_path, "sa_ssm_data.rds"))

# sa_classic_data = sensitivityAnalysis_data(model = classic, dbh0 = dbh0_ssm, sd_dbh = stanData_classic$sd_dbh, n_param = 1e3, N = 2^14,
# 	env0 = c(precip = clim_dt_ssm["precip", get(qt_opt["precip"])], tas = clim_dt_ssm["tas", get(qt_opt["tas"])], ph = 0, basalArea = 0), clim_dt = clim_dt_classic, seed = 123)

# saveRDS(sa_classic_data, paste0(tree_path, "sa_classic_data.rds"))


# setkey(sa_ssm_data, parameters, sensitivity)
# setkey(sa_classic_data, parameters, sensitivity)

# sa_ssm_data[sensitivity == "Si", sum(original), by = run][, range(V1)]
# sa_classic_data[sensitivity == "Si", sum(original), by = run][, range(V1)]


# ####! END CRASH TEST ZONE II
