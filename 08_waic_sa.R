
#### Aim of prog: Compare the SSM approach to the classic approach using the PSIS-LOO CV approach and sensitivity analysis
## Comments
# For this program, I use radial increment based on wood cores from Martínez-Sancho 2020: https://doi.org/10.1038/s41597-019-0340-y
# 	Dendroecological Collection, tree-ring and wood density data from seven tree species across Europe.
# These data have not been used to parametrise the growth model because they do not match the error structure used in our model
# these data are only available for the 7 following species:
# 	- Betula pendula
# 	- Fagus sylvatica
# 	- Picea abies
# 	- Pinus pinaster
# 	- Pinus sylvestris
# 	- Quercus petraea
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
#
## The theory of WAIC and PSIS-LOO CV is described in Vehtari.2017:
#	Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC
#	DOI: 10.1007/s11222-016-9696-4
#	https://doi.org/10.1007/s11222-016-9696-4
#
#	Note that I stupidly reprogrammed the function waic in an older version and found exactly the same results! This version uses the package loo
#
## The PSIS-LOO's abnormaly high values:
#	For some individual time series in Picea abies and Pinus sylvestris, I noted abnormaly high Pareto values. Here are, I think, the reasons:
#		- Picea abies: The concerned trees have a dbh larger than 1300 mm, while in the data used in the model, 99% of the dbh are below 647 mm
#		- Pinus sylvestris: The concerned trees all experienced extreme precip (> 2500 mm/year), while in the data 95% is below 1200 mm/year
#	This warnings are therefore not that important. Note that they also concern similar trees, in the sense that all the concerned trees for SSM
#		are also concerned trees for Classic (the opposite not being true)
#
## The theory of sensitivity analysis is described simply in Puy.2021:
#	sensobol: an R package to compute variance-based sensitivity indices
#	DOI: 10.48550/arxiv.2101.10103
#	https://arxiv.org/abs/2101.10103
#
#	Note that for a better understanding, I recommend the book of Saltelli 2008:
#	Saltelli, A.; Ratto, M.; Andres, T.; Campolongo, F.; Cariboni, J.; Gatelli, D.; Saisana, M. & Tarantola, S.
#	Global Sensitivity Analysis: The Primer
#	DOI: 10.1002/9780470725184
#	https://doi.org/10.1002/9780470725184
#
#	Note that I changed the package sensobol to the package sensitivity. This is because this package performs better with small indices
#	and is also faster. However, it is less automatised and the user manual is not very well explained... The only thing you need to know
#	is that the matrices X1 and X2 corresponds to the matrices A and B of 
#
## Climate range for SA:
#	I want to use the same climate range for both SSM and Classic SA. Indeed, it would make the comparison impossible otherwise!
#	The most logical is to use the range from the annual climate (i.e., SSM)
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500, warn = 1) # The option warn = 1 is to make the warning appears directly rahter than after the function returns.

library(sensitivity)
library(data.table)
library(tikzDevice)
library(cmdstanr)
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

## Wrap of the growth_fct_meanlog adapted to a matrix format
growth_fct_meanlog_mat = function(X, params_vec, sd_dbh, standardised_variables)
{
	varnames = c("dbh", "pr", "tas", "ph", "ba")
	if (!all(varnames %in% colnames(X)))
		stop("The colnames of the matrix X mismatch the required variables of the model")
	
	return (growth_fct_meanlog(dbh = X[, "dbh"], pr = X[, "pr"], tas = X[, "tas"], ph = X[, "ph"], basalArea = X[, "ba"],
		params = params_vec, sd_dbh = sd_dbh, standardised_variables = standardised_variables))
}

## Function to compute sensitivity of growth with respect to uncertainty in the data only
sensitivityAnalysis_data = function(model, dbh_lim, sd_dbh, clim_dt, n_param, lim_inf, lim_sup, env0 = NULL, N = 2^14, seed = NULL)
{
	## Check dbh_lim and env0
	if (!is.null(env0))
	{
		if ((env0["pr"] < clim_dt["q05", pr]) || (env0["tas"] < clim_dt["q05", tas]) || (env0["ph"] < clim_dt["q05", ph]))
			warning("Are you sure that env0 is scaled? There are low values below the 5 quantile")

		if ((env0["pr"] > clim_dt["q95", pr]) || (env0["tas"] > clim_dt["q95", tas]) || (env0["ph"] > clim_dt["q95", ph]))
			warning("Are you sure that env0 is scaled? There are large values beyond the 95 quantile")
	}

	if (length(dbh_lim) < 2)
		stop("dbh_lim must be a vector of 2 values")

	if (length(dbh_lim) > 2)
	{
		warning("Only the first two values of dbh_lim are used")
		dbh_lim = dbh_lim[1:2]
	}

	if (any(abs(dbh_lim) > 5))
		warning(paste0("Is dbh_lim scaled? dbh_lim = (", round(dbh_lim[1], 2), ", ", round(dbh_lim[2], 2), ")"))

	if (length(dbh_lim) == 2)
	{
		if (dbh_lim[1] > dbh_lim[2])
		{
			dbh_lim[1] = dbh_lim[2] - dbh_lim[1] + (dbh_lim[2] = dbh_lim[1])
			warning("The values for dbh_lim were swaped")
		}
	}
	
	## Prepare the sampling matrices
	# Common variables
	params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2",
		"competition_slope")

	if (!is.null(seed))
		set.seed(seed)
	
	if (n_param > 1)
	{
		params_values = getParams(model_cmdstan = model, params_names = params, type = "all")
	} else {
		params_values = getParams(model_cmdstan = model, params_names = params, type = "median")
		params_values = array(data = params_values, dim = c(1, 1, length(params_values)), dimnames = list("iteration", "chain", params))
		params_values = posterior::as_draws_array(params_values)
	}
	
	n_iter = model$metadata()$iter_sampling
	n_chains = model$num_chains()
	n_tot = n_iter*n_chains
	
	if (n_param > 1)
	{
		if (n_param > n_tot)
		{
			warning(paste("n_param is larger than the amount of available samples. Value set to the maximum available:", n_tot))
			n_param = n_tot
		}
		sample_ind = sample(x = 1:n_tot, size = n_param, replace = FALSE)
	} else {
		sample_ind = 1
	}

	explanatory_vars = c("dbh", "pr", "tas", "ph", "ba")
	k = length(explanatory_vars)

	## Create matrices
	sobol_mat = matrix(data = runif(n = N*2*k), nrow = N, ncol = 2*k)

	X1 = sobol_mat[, 1:k] # Correspond to matrix A in Saltelli 2008, p. 165
	colnames(X1) = explanatory_vars

	X2 = sobol_mat[, (k + 1):(2*k)] # Correspond to matrix B in Saltelli 2008, p. 165
	colnames(X2) = explanatory_vars

	# Rescale the Sobol matrix
	# --- Explanatory variables
	for (currentVar in explanatory_vars)
	{
		if (currentVar == "dbh")
		{
			X1[, "dbh"] = qunif(p = X1[, "dbh"], min = dbh_lim[1], max = dbh_lim[2])
			X2[, "dbh"] = qunif(p = X2[, "dbh"], min = dbh_lim[1], max = dbh_lim[2])
		} else {
			if (!is.null(env0))
			{
				X1[, currentVar] = qnorm(p = X1[, currentVar], mean = env0[currentVar], sd = clim_dt[currentVar, sd_sa])
				X2[, currentVar] = qnorm(p = X2[, currentVar], mean = env0[currentVar], sd = clim_dt[currentVar, sd_sa])
			} else {
				X1[, currentVar] = qunif(p = X1[, currentVar], min = clim_dt[lim_inf, get(currentVar)],
					max = clim_dt[lim_sup, get(currentVar)])
				X2[, currentVar] = qunif(p = X2[, currentVar], min = clim_dt[lim_inf, get(currentVar)],
					max = clim_dt[lim_sup, get(currentVar)])
			}
		}
	}
	
	if (any(is.na(X1)) || any(is.na(X2)))
		stop(paste("Number of NAs in Sobol's matrice:", sum(is.na(X1)) + sum(is.na(X2))))

	## Run sensitivity analysis
	ind = vector(mode = "list", length = n_param)
	
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

		# --- Sobol indices
		output_sa = sobolmartinez(model = growth_fct_meanlog_mat, X1 = X1, X2 = X2, nboot = 0, conf = 0.95,
			params_vec = params_vec, sd_dbh = sd_dbh, standardised_variables = TRUE)
		ind[[j]] = setDT(output_sa$S, keep.rownames = TRUE)

		if (any(ind[[j]][, "original"] < 0))
			warning(paste("There are negatives Si, with min value:", round(min(ind[[j]][, "original"], 7)), "for index j =", j))

		if (j %% freq_print == 0)
			print(paste0(round(100*j/n_param), "% done"))
	}
	print("100% done")

	ind = rbindlist(l = ind, idcol = "run")
	if (any(ind[, original] < 0))
		warning(paste("There are negatives Si, with min value:", round(ind[, min(original)], 7)))

	setnames(ind, old = "rn", new = "parameters")

	return(list(sa = ind, n_param = n_param))
}

#? ----------------------------------------------------------------------------------------
#* ----------------------    PART I: Compute PSIS-LOO CV and waic    ----------------------
#? ----------------------------------------------------------------------------------------
#### Load results
## Common variables
# args = c("Betula pendula", "1", "q05q95")											# NO, classic better
# args = c("Fagus sylvatica", "1", "q05q95")										# OK, ssm better
# args = c("Picea abies", "1", "q05q95")											# OK, ssm better
# args = c("Pinus pinaster", "1", "q05q95")											# OK, ssm better
# args = c("Pinus sylvestris", "1", "q05q95")										# OK, ssm better
# args = c("Quercus petraea" , "1", "q05q95")										# OK, ssm better
args = commandArgs(trailingOnly = TRUE)
for (i in seq_along(args))
	print(paste0("Arg ", i, ": <", args[i], ">"))

args[1] = paste(args[1], args[2]) #! Because of space in species' scientific name
args = args[-2]

if (length(args) != 3)
	stop("Supply in this order the species, the run_id, and the option for the sensitivity analysis", call. = FALSE)

(species = as.character(args[1]))
(run = as.integer(args[2]))
(sa_opt = as.character(args[3]))

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
scaling_dt[variable == "standBasalArea_interp", variable := "ba"] # Change the names for ease of usage

setkey(scaling_dt, type, variable)

## Create stan data to compute waic
# SSM
stanData_ssm = readRDS(paste0(tree_path, run, "_stanData.rds"))
stanData_ssm$n_data_rw = dendro[, .N]

stanData_ssm$precip_rw = dendro[, (pr - scaling_dt[.("ssm", "pr"), mu])/scaling_dt[.("ssm", "pr"), sd]]
stanData_ssm$tas_rw = dendro[, (tas - scaling_dt[.("ssm", "tas"), mu])/scaling_dt[.("ssm", "tas"), sd]]
stanData_ssm$ph_rw = dendro[, (ph - scaling_dt[.("ssm", "ph"), mu])/scaling_dt[.("ssm", "ph"), sd]]
stanData_ssm$standBasalArea_rw = dendro[, standBasalArea]

stanData_ssm$dbh_rw = dendro[, current_dbh/stanData_ssm$sd_dbh]
stanData_ssm$ring_width = dendro[, dbh_increment_in_mm/stanData_ssm$sd_dbh]

# Classic
stanData_classic = readRDS(paste0(tree_path, run, "_stanData_classic.rds"))
stanData_classic$n_data_rw = dendro[, .N]

stanData_classic$precip_rw = dendro[, (pr - scaling_dt[.("classic", "pr"), mu])/scaling_dt[.("classic", "pr"), sd]]
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
saveRDS(loo_ssm, paste0(tree_path, "loo_ssm.rds"))

r_eff_classic = loo::relative_eff(loglik_classic$draws("log_lik"), cores = 8)
loo_classic = loo::loo(x = loglik_classic$draws("log_lik"), r_eff = r_eff_classic, cores = 8)
waic_classic = loo::waic(x = loglik_classic$draws("log_lik"), cores = 8)
saveRDS(loo_classic, paste0(tree_path, "loo_classic.rds"))

loo::loo_compare(list(loo_ssm, loo_classic))
loo::loo_compare(list(waic_ssm, waic_classic))

#? -----------------------------------------------------------------------------------------
#* ----------------------    PART II: Compute sensitivity analysis    ----------------------
#? -----------------------------------------------------------------------------------------
#### Climate range (see comments at the beginning)
if ((file.exists("./speciesInformations.rds")) && (file.exists("./speciesInformations_runs.rds")))
{
	info = readRDS("./speciesInformations.rds")
} else {
	ls_info = infoSpecies()
}

info = data.table::melt(data = info[species], measure = patterns("^dbh_", "^tas", "^pr_", "^ph_", "^ba_"),
	value.name = c("dbh", "tas", "pr", "ph", "ba"))

info[, variable := NULL]
info[, qtile := c("min", "q025", "q05", "med", "q95", "q975", "max")]

for (i in 2:info[, .N])
{
	previousLine = info[i - 1, c(dbh, tas, pr, ph, ba)]
	currentLine = info[i, c(dbh, tas, pr, ph, ba)]

	if (any(currentLine - previousLine < 0))
		stop("The qtile are wrongly assigned, check clim_dt")
}

setkey(info, qtile)

info_ssm = info[, .(qtile, dbh, tas, pr, ph, ba)]
info_classic = info[, .(qtile, dbh, tas, pr, ph, ba)]

for (currentVar in c("tas", "pr", "ph", "ba")) # Same data table, but scaled differently!
{
	info_ssm[, (currentVar) := (get(currentVar) - scaling_dt[.("ssm", currentVar), mu])/scaling_dt[.("ssm", currentVar), sd]]
	info_classic[, (currentVar) := (get(currentVar) - scaling_dt[.("classic", currentVar), mu])/scaling_dt[.("classic", currentVar), sd]]
}

#### Sensitivity of growth with respect to data only
## Common variables (note that dbh_lim and sd_dbh are the same between SSM and Classic)
if (sa_opt == "q025q975")
{
	lim_inf = "q025"
	lim_sup = "q975"
} else if (sa_opt == "min_max")
{
	lim_inf = "min"
	lim_sup = "max"
} else {
	lim_inf = "q05"
	lim_sup = "q95"
	if (sa_opt != "q05q95")
		warning("lim_inf and lim_sup are set to default values: q05 and q95, respectively")
}

dbh_lim = c(info_ssm[c(lim_inf, lim_sup), dbh])/stanData_ssm$sd_dbh
sd_dbh = stanData_ssm$sd_dbh
n_param = 500

## SSM
sa_ssm_data = sensitivityAnalysis_data(model = ssm, dbh_lim = dbh_lim, sd_dbh = sd_dbh, n_param = n_param, lim_inf = lim_inf,
	lim_sup = lim_sup, N = 2^19, clim_dt = info_ssm, seed = 123)

saveRDS(sa_ssm_data, paste0(tree_path, "sa_ssm_data_", lim_inf, "_", lim_sup, "_nParams=", n_param, ".rds"))

## Classic
sa_classic_data = sensitivityAnalysis_data(model = classic, dbh_lim = dbh_lim, sd_dbh = sd_dbh, n_param = n_param, lim_inf = lim_inf,
	lim_sup = lim_sup, N = 2^19, clim_dt = info_classic, seed = 123)

saveRDS(sa_classic_data, paste0(tree_path, "sa_classic_data_", lim_inf, "_", lim_sup, "_nParams=", n_param, ".rds"))

print(sa_ssm_data$sa[, sum(original), by = run][, range(V1)])
print(sa_classic_data$sa[, sum(original), by = run][, range(V1)])

####! CRASH TEST ZONE
f = function(dbh, pr, tas, ph, params)
{
	results = params["averageGrowth"] +
		params["dbh_slope"] * dbh + params["dbh_slope2"] * dbh^2 +
		params["pr_slope"] * pr + params["pr_slope2"] * pr^2 +
		params["tas_slope"] * tas + params["tas_slope2"] * tas^2 +
		params["ph_slope"] * ph + params["ph_slope2"] * ph^2
	return(exp(results))
}

f_transform = function(dbh, pr, tas, ph, params)
{
	results = params[["normalising"]] * exp(toTry[["intercept"]]) *
		ifelse(params[["original"]]["dbh_slope2"] < 0, dnorm(dbh, params[["mu"]]["dbh"], params[["sigma"]]["dbh"]),
			exp(toTry[["original"]]["dbh_slope"]*dbh + toTry[["original"]]["dbh_slope2"]*dbh^2)) *
		ifelse(params[["original"]]["pr_slope2"] < 0, dnorm(pr, params[["mu"]]["pr"], params[["sigma"]]["pr"]),
			exp(toTry[["original"]]["pr_slope"]*pr + toTry[["original"]]["pr_slope2"]*pr^2)) *
		ifelse(params[["original"]]["tas_slope2"] < 0, dnorm(tas, params[["mu"]]["tas"], params[["sigma"]]["tas"]),
			exp(toTry[["original"]]["tas_slope"]*tas + toTry[["original"]]["tas_slope2"]*tas^2)) *
		ifelse(params[["original"]]["ph_slope2"] < 0, dnorm(ph, params[["mu"]]["ph"], params[["sigma"]]["ph"]),
			exp(toTry[["original"]]["ph_slope"]*ph + toTry[["original"]]["ph_slope2"]*ph^2))
	return(results)
}

testParams = c(averageGrowth = 2.3,
	dbh_slope = 1.4, dbh_slope2 = 0.25,
	pr_slope = 0.24, pr_slope2 = -1.2,
	tas_slope = 0.8, tas_slope2 = 0.05,
	ph_slope = -1.23, ph_slope2 = 0.87)

toTry = gaussStyle(params = testParams)

f(0.8, 1.2, -0.87, 0.56, testParams)
f_transform(0.8, 1.2, -0.87, 0.56, toTry)
