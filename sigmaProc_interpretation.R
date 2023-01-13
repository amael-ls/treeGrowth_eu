
#### Aim of prog: Comparing sigmaProc (SSM approach) vs sigmaProc (classic approach)
# Comments:
# The package lognorm is based on http://www.m-hikari.com/ams/ams-2013/ams-125-128-2013/39511.html
#	"WKB Approximation for the Sum of Two Correlated Lognormal Random Variables".
#	I use this package to approximate the sum of lognormal random variables by a single lognormal distribution. The article looks fine,
#	although the author did not check the relative error of his method, only the absolute error. The relative error shows that his method
#	is unable to correclty fit the tail or the head of the distribution. Other methods have been developed where the precision in the tail
#	or head is also discussed in https://ieeexplore.ieee.org/document/1578407 "A Flexible Lognormal Sum Approximation Method".
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(lognorm)

#### Tool functions
meanlog_fct = function(dbh0, precip, temperature, ph, standBasalArea, averageGrowth, dbh_slope, dbh_slope2, pr_slope, pr_slope2,
	tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope)
{
	return(averageGrowth + dbh_slope*dbh0 + dbh_slope2*dbh0^2 + pr_slope*precip + pr_slope2*precip^2 +
		tas_slope*temperature + tas_slope2*temperature^2 + ph_slope*ph + ph_slope2*ph^2 + competition_slope*standBasalArea)
}

source("./toolFunctions.R")

#### Common variables
## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

n_runs = 4 # Number of runs used in growth_subsample.R
threshold_indiv = 12000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2023/01/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed_ssm := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_de-fr-sw_12000_main.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
threshold_time = as.Date("2022/12/14") # For Fagus sylvatica
infoSpecies[, processed_classic := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_main_classic.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed_ssm) & (processed_classic)]



#### For loop on processed species
for (species in infoSpecies[, speciesName_sci])
{
	#TODO Run a function that will compute the approximation of the sum of lognormal RVs and then compare it to the classical approach
}









########! WAINTING FOR ALL THE RESULTS, PREPARING THE PROGRAM (WITH TEMPORARY RESULTS, BUT THE COLUMNS/STRUCTURES DO NOT CHANGE)
#### Load data
## General data
# Trees
treeData_classic = readRDS("./Fagus sylvatica/1_treeData_classic.rds")
treeData_ssm = readRDS("./Fagus sylvatica/1_treeData.rds") # Should be the same as classic data (currently not because 8000 vs 12000 individuals)
treeData = readRDS("./Fagus sylvatica/1_treeData.rds") #! This should be the only treeData loaded as ssm and classic data are the same!

# Climate
clim_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))
climate = readRDS(paste0(clim_folder, "europe_reshaped_climate.rds"))

# Indices
indices = readRDS("./Fagus sylvatica/1_indices.rds")[["indices_avg"]]
indices[, year := NULL]
nClim = indices[, sum(year_end - year_start + 1)]

dataClim = setnames(data.table(matrix(data = 0, nrow = nClim, ncol = 6)), names(climate))
dataClim[, plot_id := as.character(plot_id)]
cols = names(climate)
start = 1

for (i in seq_len(indices[, .N]))
{
	end = start + indices[i , year_end - year_start]
	dataClim[start:end, c(cols) := climate[indices[i, index_clim_start]:indices[i, index_clim_end]]]
	start = end + 1

	if (i %% 100 == 0)
		print(paste0(i*100/indices[, .N], "% done"))
}
print("100% done")

setkey(dataClim, plot_id, year)

# Scalings
scaling_dbh = readRDS("Fagus sylvatica/1_dbh_normalisation.rds")
setkey(scaling_dbh, variable)
scaling_ba = readRDS("Fagus sylvatica/1_ba_normalisation.rds")
setkey(scaling_ba, variable)
scaling_clim = readRDS("Fagus sylvatica//1_climate_normalisation.rds")
setkey(scaling_clim, variable)

## Classic approach
classic = readRDS("./Fagus sylvatica/growth-run=1-2022-12-15_11h06_de-fr-sw_8000_main_classic.rds")
# classic$print(c("lp__", "averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
# 	"ph_slope", "ph_slope2", "competition_slope", "etaObs", "proba", "sigmaProc"), max_rows = 20)

## SSM
ssm = readRDS("./Fagus sylvatica/growth-run=1-2023-01-12_02h02_de-fr-sw_12000_main.rds")
# ssm$print(c("lp__", "averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
# 	"ph_slope", "ph_slope2", "competition_slope", "etaObs", "proba", "sigmaProc"), max_rows = 20)

#### Compute the approximation of the lognormal 
## Common variables
selected_plot_id = "france_501878"
selected_tree_id = "14"
selected_year = 2010

treeData[.(selected_plot_id, selected_tree_id, selected_year)]

data[["scaling"]]["sd_dbh_orig"]

init_dbh = unlist(treeData[.(selected_plot_id, selected_tree_id, selected_year), dbh])/scaling_dbh["dbh", sd]
current_dbh = init_dbh



#! RESTART HERE ----------------------------------------------------------------------------------------------------------------
climCols = names(dataClim)[!(names(dataClim) %in% c("plot_id", "year"))]


dataClim[.(selected_plot_id, selected_year), ..climCols]

temperatures = (unlist(treeData[.(selected_plot_id, selected_tree_id), ..temperature_cols]) - data[["scaling"]]["mu_temp"])/
	data[["scaling"]]["sd_temp"]






nGrowth = data[["infos"]]["n_annual_growth_per_indiv"]
nSim = 1e6

meanlog = numeric(nGrowth)
dt_sample = setNames(data.table(matrix(data = 0, nrow = nSim, ncol = nGrowth)), paste0("sample", 1:nGrowth))

## Parameters value (averaged estimates)
beta0 = mean(ssm$draws("beta0"))
beta1 = mean(ssm$draws("beta1"))
beta2 = mean(ssm$draws("beta2"))
beta3 = mean(ssm$draws("beta3"))
beta4 = mean(ssm$draws("beta4"))

sdlog = rep(mean(ssm$draws("sigmaProc")), nGrowth)

for (i in 1:nGrowth)
{
	dbhCol = paste0("dbh", i)
	meanlog[i] = meanlog_fct(dbh0 = current_dbh, temperature = temperatures[i], beta0, beta1, beta2, beta3, beta4)
	sampleCol = paste0("sample", i)
	dt_sample[, c(sampleCol) := rlnorm(n = nSim, meanlog = meanlog[i], sdlog = sdlog[i])]
	current_dbh = current_dbh + unlist(dt_sample[, lapply(.SD, mean), .SDcols = c(sampleCol)])
}

## Approximated lognormal distribution from the sum of lognormal random variables
(coefSum = estimateSumLognormal(mu = meanlog, sigma = sdlog))

## Comparison with the parameters from the classic approach
mean(classic$draws("sigmaProc"))
avg_temperature = (unlist(treeData[.(selected_plot_id, selected_tree_id), avg_temperature1]) - data[["scaling"]]["mu_temp"])/
	data[["scaling"]]["sd_temp"]

meanlog_fct(dbh0 = init_dbh, temperature = avg_temperature, beta0 = mean(classic$draws("beta0")), beta1 = mean(classic$draws("beta1")),
	beta2 = mean(classic$draws("beta2")), beta3 = mean(classic$draws("beta3")), beta4 = mean(classic$draws("beta4"))) +
	log(data[["infos"]]["delta_t"])
	# The + log(delta_t) is there because it is dbh(t + delta_t) = dbh(t) + delta_t * logN(...) in the classic approach 
