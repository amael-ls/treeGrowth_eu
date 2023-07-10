
#! A COPY OF compareModels.R with modifications to use waic.stan
#### Aim of prog: Compare the SSM approach to the classic approach using the WAIC
## Comments
# The WAIC is from Vehtari.2017:
#	Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC
#	DOI: 10.1007/s11222-016-9696-4
#	https://doi.org/10.1007/s11222-016-9696-4


#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

# library(microbenchmark)
library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)
library(Rmpfr)
library(terra)
library(qgam)

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

## Function to compute log predictive density
waic = function(loglik_ssm, loglik_classic, n_chains, n_iter)
{
	draws_ssm = loglik_ssm$draws()
	draws_classic = loglik_classic$draws()

	S = n_chains*n_iter
	
	mean_draws_ssm = apply(X = exp(draws_ssm), MARGIN = 3, FUN = mean) # The exp is necessary because stan takes the log of proba distrib fct
	mean_draws_classic = apply(X = exp(draws_classic), MARGIN = 3, FUN = mean)

	if ((min(mean_draws_ssm) == 0) || (min(mean_draws_classic) == 0))
		print("Zeros are present but not removed, expect a -Inf!")
	
	lpd_hat_ssm = sum(log(mean_draws_ssm))
	lpd_hat_classic = sum(log(mean_draws_classic))

	n_indiv_rw = length(mean_draws_ssm)
	p_hat_ssm = numeric(n_indiv_rw)
	p_hat_classic = numeric(n_indiv_rw)

	for (measure in 1:n_indiv_rw)
	{
		temporary = (draws_ssm[, , measure] - 1/S*sum(draws_ssm[, , measure]))^2
		p_hat_ssm[measure] = 1/(S - 1) * sum(temporary)

		temporary = (draws_classic[, , measure] - 1/S*sum(draws_classic[, , measure]))^2
		p_hat_classic[measure] = 1/(S - 1) * sum(temporary)
	}
	print(paste("min_ssm = ", min(p_hat_ssm)))
	print(paste("max_ssm = ", max(p_hat_ssm)))
	print(paste("min_classic = ", min(p_hat_classic)))
	print(paste("max_classic = ", max(p_hat_classic)))

	return(c(waic_ssm = lpd_hat_ssm - sum(p_hat_ssm), waic_classic = lpd_hat_classic - sum(p_hat_classic)))
}


#? --------------------------------------------------------------------------------------------------
#* ----------------------    PART I: Compute quantiles GAM on the residuals    ----------------------
#? --------------------------------------------------------------------------------------------------
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
dbh_scaling_ssm = readRDS(paste0(tree_path, run, "_dbh_normalisation.rds"))
ba_scaling_ssm = readRDS(paste0(tree_path, run, "_ba_normalisation.rds"))
climate_scaling_ssm = readRDS(paste0(tree_path, run, "_climate_normalisation.rds"))
ph_scaling_ssm = readRDS(paste0(tree_path, run, "_ph_normalisation.rds"))

dbh_scaling_classic = readRDS(paste0(tree_path, run, "_dbh_normalisation_classic.rds"))
ba_scaling_classic = readRDS(paste0(tree_path, run, "_ba_normalisation_classic.rds"))
climate_scaling_classic = readRDS(paste0(tree_path, run, "_climate_normalisation_classic.rds"))
ph_scaling_classic = readRDS(paste0(tree_path, run, "_ph_normalisation_classic.rds"))

scaling_ssm = rbindlist(list(dbh_scaling_ssm, ba_scaling_ssm, climate_scaling_ssm, ph_scaling_ssm))
scaling_classic = rbindlist(list(dbh_scaling_classic, ba_scaling_classic, climate_scaling_classic, ph_scaling_classic))

scaling_dt = rbindlist(list(ssm = scaling_ssm, classic = scaling_classic), idcol = "type")
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

stanData_classic$precip_rw = dendro[, (pr - scaling_dt[.("classic", "pr_avg"), mu])/scaling_dt[.("classic", "pr_avg"), sd]]
stanData_classic$tas_rw = dendro[, (tas - scaling_dt[.("classic", "tas_avg"), mu])/scaling_dt[.("classic", "tas_avg"), sd]]
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

waic_hat = waic(loglik_ssm = loglik_ssm, loglik_classic = loglik_classic, n_chains = n_chains, n_iter = n_iter)
