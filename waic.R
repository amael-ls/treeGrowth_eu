
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
#	Note that for a better understanding, the book Sensitivity Analysis: The Primer (Saltelli, 2008) is a good option.
#	"Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index" xz(Saltell.2010) gives an overview 
#	of Total Sensitivity Indices.
#	DOI: 10.1016/j.cpc.2009.09.018
#	https://www.sciencedirect.com/science/article/abs/pii/S0010465509003087?via%3Dihub
#
## Explanations on the tree rings data
#	The data contains for each year the dbh increment, in mm, for that year. In the following example, the dbh increment from the 1st January
#		1980 to 31st December 1980 is 0.16 mm.
# 
#		   plot_id  tree_id longitude latitude [...] dbh  year dbh_increment_in_mm starting_dbh       ph      pr      tas
#		1:  ATFS13 ATFS1304  14.06788 46.49724 [...] 200  1980                0.16       154.64 5.311581 2063.37 4.000000
#		2:  ATFS13 ATFS1304  14.06788 46.49724 [...] 200  1981                0.38       154.64 5.311581 1692.40 4.775000
#		[...]
#		36: ATFS13 ATFS1304  14.06788 46.49724 [...] 200  2015                3.14       154.64 5.311581 2030.09 7.375000
#	
#		and the increment in 2015 is 3.14 mm. The column dbh corresponds to the diameter measured one year after the last,
#		i.e., in 2016 in this example.
#
#	I defined the starting_dbh as the dbh the tree would have at the first year of the record (here, 1980). I computed it by substracting
#		the sum of increments to the dbh. Note that it neglects the growth that occured the year of dbh measurement (here, 2016) but it
#		should not affect the results: Growth is a continuous function of dbh, and its derivative with respect to dbh (which is dbh_slope)
#		is 'small' (whatever that means!)

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
# args = c("Betula pendula", "1")											# NO, classic better
# args = c("Fagus sylvatica", "1")											# OK, ssm better
# args = c("Picea abies", "1")												# OK, ssm better
# args = c("Pinus pinaster", "1")											# OK, ssm better
# args = c("Pinus sylvestris", "1")											# OK, ssm better
# args = c("Quercus petraea" , "1")											# OK, ssm better
args = commandArgs(trailingOnly = TRUE)
for (i in seq_along(args))
	print(paste0("Arg ", i, ": <", args[i], ">"))

args[1] = paste(args[1], args[2])
args = args[-2]

if (length(args) != 2)
	stop("Supply in this order the species and the run_id", call. = FALSE)

(species = as.character(args[1]))
(run = as.integer(args[2]))

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

#### Sensitivity of growth with respect to data only
## SSM
sa_ssm_data = sensitivityAnalysis_data(model = ssm, dbh0 = quantile(stanData_ssm$dbh_init, c(0.05, 0.95))/stanData_ssm$sd_dbh,
	sd_dbh = stanData_ssm$sd_dbh, n_param = 500, N = 2^16, clim_dt = clim_dt_ssm, seed = 123)

saveRDS(sa_ssm_data, paste0(tree_path, "sa_ssm_data.rds"))

## Classic
sa_classic_data = sensitivityAnalysis_data(model = classic, dbh0 = quantile(stanData_ssm$dbh_init, c(0.05, 0.95))/stanData_ssm$sd_dbh,
	sd_dbh = stanData_classic$sd_dbh, n_param = 500, N = 2^16, clim_dt = clim_dt_classic, seed = 123)

saveRDS(sa_classic_data, paste0(tree_path, "sa_classic_data.rds"))


print(sa_ssm_data$sa[sensitivity == "Si", sum(original), by = run][, range(V1)])
print(sa_classic_data$sa[sensitivity == "Si", sum(original), by = run][, range(V1)])
