
#### Aim of prog: Comparing sigmaProc (SSM approach) vs sigmaProc (classic approach)
## Comments:
# The package lognorm is based on http://www.m-hikari.com/ams/ams-2013/ams-125-128-2013/39511.html
#	"WKB Approximation for the Sum of Two Correlated Lognormal Random Variables".
#	I use this package to approximate the sum of lognormal random variables by a single lognormal distribution. The article looks fine,
#	although the author did not check the relative error of his method, only the absolute error. The relative error shows that his method
#	is unable to correclty fit the tail or the head of the distribution. Other methods have been developed where the precision in the tail
#	or head is also discussed in https://ieeexplore.ieee.org/document/1578407 "A Flexible Lognormal Sum Approximation Method".
#
# How do I compare the two approaches:
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










#### Load data
## General data
# Trees
treeData = readRDS("./Fagus sylvatica/1_treeData.rds")

# Climate
clim_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

climate = readRDS(paste0(clim_folder, "europe_reshaped_climate.rds"))

# Basal area (interpolated data)
standBasalArea_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(standBasalArea_folder))
	stop(paste0("Folder\n\t", standBasalArea_folder, "\ndoes not exist"))

standBasalArea = readRDS(paste0(standBasalArea_folder, "europe_reshaped_standBasalArea.rds"))

# Soil data
soil_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(soil_folder))
	stop(paste0("Folder\n\t", soil_folder, "\ndoes not exist"))

soil = readRDS(paste0(soil_folder, "europe_reshaped_soil.rds"))
setkey(soil, plot_id)

# Merge climate  stand basal area data
climate = climate[standBasalArea, on = c("plot_id", "year")]

# Indices
indices = readRDS("./Fagus sylvatica/1_indices.rds")[["indices_avg"]]
indices[, year := NULL]
nClim = indices[, sum(year_end - year_start + 1)]

dataClim = setnames(data.table(matrix(data = 0, nrow = nClim, ncol = ncol(climate))), names(climate))
dataClim[, plot_id := as.character(plot_id)]
cols = names(climate)
start = 1

for (i in seq_len(indices[, .N]))
{
	end = start + indices[i , year_end - year_start]
	dataClim[start:end, c(cols) := climate[indices[i, index_clim_start]:indices[i, index_clim_end]]]
	start = end + 1

	if (i %% 100 == 0)
		print(paste0(round(i*100/indices[, .N], 2), "% done"))
}
print("100% done")

setkey(dataClim, plot_id, year)

# Scalings
scaling_dbh = readRDS("Fagus sylvatica/1_dbh_normalisation.rds")
setkey(scaling_dbh, variable)
scaling_ba = readRDS("Fagus sylvatica/1_ba_normalisation.rds")
setkey(scaling_ba, variable)
scaling_clim = readRDS("Fagus sylvatica/1_climate_normalisation.rds")
setkey(scaling_clim, variable)
scaling_ph = readRDS("Fagus sylvatica/1_ph_normalisation.rds")
setkey(scaling_ph, variable)

## Classic approach
classic = readRDS("./Fagus sylvatica/growth-run=1-2023-01-20_20h22_de-fr-sw_12000_main_classic.rds")

## SSM
ssm = readRDS("./Fagus sylvatica/growth-run=1-2023-01-12_02h02_de-fr-sw_12000_main.rds")

#### Compute the approximation of the lognormal 
## Common variables
# selected_plot_id = "france_501878"
# selected_tree_id = "14"
# selected_year = 2010

selected_plot_id = "germany_12591_3"
selected_tree_id = "10"
selected_year = 2000

treeData[.(selected_plot_id, selected_tree_id, selected_year)]
treeData[.(selected_plot_id, selected_tree_id)]

init_dbh = unlist(treeData[.(selected_plot_id, selected_tree_id, selected_year), dbh])/scaling_dbh["dbh", sd]
current_dbh = init_dbh

nGrowth = treeData[.(selected_plot_id, selected_tree_id, selected_year), deltaYear]
selected_years_clim = selected_year:(selected_year + nGrowth)

temperatures = (dataClim[.(selected_plot_id, selected_years_clim), tas] - scaling_clim["tas", mu])/scaling_clim["tas", sd]
precipitations = (dataClim[.(selected_plot_id, selected_years_clim), pr] - scaling_clim["pr", mu])/scaling_clim["pr", sd]
basalAreas = (dataClim[.(selected_plot_id, selected_years_clim), standBasalArea_interp] - scaling_ba["standBasalArea_interp", mu])/
	scaling_ba["standBasalArea_interp", sd]
ph = (soil[selected_plot_id, ph] - scaling_ph["ph", mu])/scaling_ph["ph", sd]
nSim = 1e6

plot(selected_years_clim, temperatures, pch = 19, ylim = c(min(c(temperatures, precipitations)), max(c(temperatures, precipitations))),
	col = "#CD212A")
points(selected_years_clim, precipitations, pch = 19, col = "#122A51")
lines(selected_years_clim, temperatures, col = "#CD212A")
lines(selected_years_clim, precipitations, col = "#122A51")
abline(h = mean(temperatures), lwd = 2, col = "#CD212A")
abline(h = mean(precipitations), lwd = 2, col = "#122A51")
meanlog = numeric(nGrowth)
dt_sample = setNames(data.table(matrix(data = 0, nrow = nSim, ncol = nGrowth)), paste0("sample", 1:nGrowth))

## Parameters value (averaged estimates)
averageGrowth = mean(ssm$draws("averageGrowth"))
dbh_slope = mean(ssm$draws("dbh_slope"))
dbh_slope2 = mean(ssm$draws("dbh_slope2"))
pr_slope = mean(ssm$draws("pr_slope"))
pr_slope2 = mean(ssm$draws("pr_slope2"))
tas_slope = mean(ssm$draws("tas_slope"))
tas_slope2 = mean(ssm$draws("tas_slope2"))
ph_slope = mean(ssm$draws("ph_slope"))
ph_slope2 = mean(ssm$draws("ph_slope2"))
competition_slope = mean(ssm$draws("competition_slope"))

sdlog = rep(mean(ssm$draws("sigmaProc")), nGrowth)

for (i in 1:nGrowth)
{
	dbhCol = paste0("dbh", i)
	meanlog[i] = meanlog_fct(dbh0 = current_dbh, precip = precipitations[i], temperature = temperatures[i], ph = ph,
		standBasalArea = basalAreas[i], averageGrowth, dbh_slope, dbh_slope2, pr_slope, pr_slope2, tas_slope, tas_slope2,
		ph_slope, ph_slope2, competition_slope)
	sampleCol = paste0("sample", i)
	dt_sample[, c(sampleCol) := rlnorm(n = nSim, meanlog = meanlog[i], sdlog = sdlog[i])]
	current_dbh = current_dbh + unlist(dt_sample[, lapply(.SD, mean), .SDcols = c(sampleCol)])
}

## Approximated lognormal distribution from the sum of lognormal random variables
(coefSum = estimateSumLognormal(mu = meanlog, sigma = sdlog))

sd_dbh = scaling_dbh["dbh", sd]
avgGrowth = exp(meanlog + sdlog^2/2)

cor(precipitations[1:nGrowth], avgGrowth)
cor(temperatures[1:nGrowth], avgGrowth)

## Comparison with the parameters from the classic approach
mean(classic$draws("sigmaProc"))
avg_temperature = mean(temperatures)
avg_precipitations = mean(precipitations)
avg_basalArea = mean(basalAreas)

meanlog_classic = meanlog_fct(dbh0 = init_dbh, precip = avg_precipitations, temperature = avg_temperature, ph = ph,
		standBasalArea = avg_basalArea, mean(classic$draws("averageGrowth")), mean(classic$draws("dbh_slope")),
		mean(classic$draws("dbh_slope2")), mean(classic$draws("pr_slope")), mean(classic$draws("pr_slope2")), mean(classic$draws("tas_slope")),
		mean(classic$draws("tas_slope2")), mean(classic$draws("ph_slope")), mean(classic$draws("ph_slope2")),
		mean(classic$draws("competition_slope")))
	# Add + log(nGrowth) if simulating over many years, dbh(t + nGrowth) = dbh(t) + nGrowth * logN(...) in the classic approach

curve(dlnorm(x, meanlog = coefSum["mu"] + log(sd_dbh) - log(nGrowth), sdlog = coefSum["sigma"]), lwd = 2, to = 30,
	xlab = "Averaged annual growth (mm/yr)", ylab = "Posterior", col = "#122A51")
curve(dlnorm(x, meanlog = meanlog_classic + log(sd_dbh), sdlog = mean(classic$draws("sigmaProc"))),
	col = "#A21121", lwd = 2, add = TRUE)
legend(x = "topright", legend = c("SSM - approx", "Classic"), col = c("#122A51", "#A21121"), lwd = 2)


exp(coefSum["mu"] + log(sd_dbh) - log(nGrowth) + coefSum["sigma"]^2/2)
exp(meanlog_classic + log(sd_dbh) + mean(classic$draws("sigmaProc"))^2/2)


aa2 = rlnorm(1e6, meanlog = meanlog[2] + log(sd_dbh), sdlog = sdlog[2])
aa3 = rlnorm(1e6, meanlog = meanlog[3] + log(sd_dbh), sdlog = sdlog[3])
aa4 = rlnorm(1e6, meanlog = meanlog[4] + log(sd_dbh), sdlog = sdlog[4])
aa5 = rlnorm(1e6, meanlog = meanlog[5] + log(sd_dbh), sdlog = sdlog[5])
aa6 = rlnorm(1e6, meanlog = meanlog[6] + log(sd_dbh), sdlog = sdlog[6])
aa7 = rlnorm(1e6, meanlog = meanlog[7] + log(sd_dbh), sdlog = sdlog[7])

ss = aa2 + aa3 + aa4 + aa5 + aa6 + aa7
ss = ss[ss < quantile(ss, seq(0, 1, 0.05))["95%"]]
hist(ss, prob = TRUE)
(aa_aprox = estimateSumLognormal(mu = meanlog[2:7], sigma = sdlog[2:7]))
curve(dlnorm(x, meanlog = aa_aprox["mu"] + log(sd_dbh), aa_aprox["sigma"]), add = TRUE)
