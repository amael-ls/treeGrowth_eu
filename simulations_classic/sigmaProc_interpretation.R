
#### Aim of prog: Comparing sigmaProc (SSM approach) vs sigmaProc (classic approach)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(lognorm)

#### Tool functions
meanlog_fct = function(dbh0, temperature, beta0, beta1, beta2, beta3, beta4)
	return(beta0 + beta1*dbh0 + beta2*dbh0^2 + beta3*temperature + beta4*temperature^2)

#### Load data
## General data (list, same for both models)
data = readRDS("./dummyData.rds")

## Classic approach
classic = readRDS("./indiv=2000_plot=800_measurements=2_deltaT=5_results.rds")

## SSM
ssm = readRDS("../simulations_aubry-kientz/indiv=2000_plot=800_measurements=2_deltaT=5_results.rds")

#### Compute the approximation of the lognormal 
## Common variables
setkey(data[["treeData"]], plot_id, tree_id)

selected_plot_id = 1
selected_tree_id = 1

init_dbh = unlist(data[["treeData"]][.(selected_plot_id, selected_tree_id), dbh1])/data[["scaling"]]["sd_dbh_orig"]
current_dbh = init_dbh
temperature_cols = colnames(data[["treeData"]])[stri_detect(colnames(data[["treeData"]]), regex = "temperature[:digit:]")]
temperature_cols = temperature_cols[!stri_detect(temperature_cols, regex = "avg_temperature[:digit:]")]
temperatures = (unlist(data[["treeData"]][.(selected_plot_id, selected_tree_id), ..temperature_cols]) - data[["scaling"]]["mu_temp"])/
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

(coefSum = estimateSumLognormal(mu = meanlog, sigma = sdlog))

mean(classic$draws("sigmaProc"))
avg_temperature = (unlist(data[["treeData"]][.(selected_plot_id, selected_tree_id), avg_temperature1]) - data[["scaling"]]["mu_temp"])/
	data[["scaling"]]["sd_temp"]

meanlog_fct(dbh0 = init_dbh, temperature = avg_temperature, beta0 = mean(classic$draws("beta0")), beta1 = mean(classic$draws("beta1")),
	beta2 = mean(classic$draws("beta2")), beta3 = mean(classic$draws("beta3")), beta4 = mean(classic$draws("beta4"))) +
	log(data[["infos"]]["delta_t"])
	# The + log(delta_t) is there because it is dbh(t + delta_t) = dbh(t) + delta_t * logN(...) in the classic approach 
