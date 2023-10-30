
#### Aim of prog: validation of my model using independent data from Spain.
## Comments:
#	The first part is dedicated to load, subsample, and format the different data (dbh, climate, soil, basal area)
#	The second part is dedicated to copmute the residuals
#
## Explanations:
#	The data are from the Spanish NFI number 2 (1996) and number 4 (2011). The projection used is EPSG:23030, also called ED50 / UTM zone 30N,
#		https://epsg.io/23030
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(terra)

#### Tool functions
source("toolFunctions.R")

## Function to compute growth, the data table must be sorted by year within tree id and plot id
computeDiametralGrowth = function(dt, col = "growth", byCols = c("plot_id", "tree_id"))
{
	if (!all(c("dbh", "year", byCols) %in% names(dt)))
		stop(paste("The data table must contains at least contains the columns dbh, year,", paste(byCols, collapse = ", ")))
	while (col %in% names(dt))
	{
		newcol = paste0(col, rnorm(1))
		warning(paste0("The name `", col, "` is already used in the data table. The result was stored in col `", newcol, "` instead"))
		col = newcol
	}
	dt[, (col) := (data.table::shift(dbh, n = 1, type = "lead", fill = NA) - dbh)/
		(data.table::shift(year, n = 1, type = "lead", fill = NA) - year), by = byCols]

	if (!("deltaYear" %in% names(dt)))
		dt[, deltaYear := data.table::shift(year, n = 1, type = "lead", fill = NA) - year, by = byCols]

	if (!("increment" %in% names(dt)))
		dt[, increment := data.table::shift(dbh, n = 1, type = "lead", fill = NA) - dbh, by = byCols]
}

######## Part I: Load and subsample data
#### Common variables
args = commandArgs(trailingOnly = TRUE)
args = c("Picea abies", "1")
if (length(args) != 2)
	stop("Supply the species_id and run_id!", call. = FALSE)

species = as.character(args[1])
run = as.integer(args[2])

tree_path = paste0("./", species, "/")
if (!dir.exists(tree_path))
	stop(paste0("Path not found for species <", species, ">."))

set.seed(1969 - 08 - 18) # Woodstock seed

#### Trees
## Load
plotInfo = fread("/bigdata/Inventories/ES IFN/Original data/Adult_Recruitment_SFI24/plot24.csv")
plotInfo[, V1 := NULL]
setnames(x = plotInfo, old = c("IDPC4", "plotcode", "CX", "CY", "year2", "year4"),
	new = c("plot_id_4", "plot_id_2", "lon", "lat", "year_2", "year_4"))

speciesInfo = fread("/bigdata/Inventories/ES IFN/Original data/Adult_Recruitment_SFI24/species.csv")
speciesInfo = unique(speciesInfo[, .(Code, Code3, Code2, Nombre)])

tree_data = fread("/bigdata/Inventories/ES IFN/Original data/Adult_Recruitment_SFI24/tree24_dbh.csv")
tree_data = tree_data[!is.na(dbh2) & !is.na(dbh4)]
tree_data = na.omit(tree_data)
tree_data = tree_data[(alive2 == "alive") & (alive4 == "alive")]
tree_data[, c("V1", "alive2", "alive4") := NULL]

setnames(x = tree_data, old = c("ID_Pma3c", "IDPC4", "Plotcode4", "sppcompa", "dbh2", "dbh4"),
	new = c("tree_id", "plot_id_4", "plot_id_2", "species_code", "dbh_2", "dbh_4"))

## Merge tree data with plot data and reshape
tree_data = plotInfo[tree_data, on = c("plot_id_4", "plot_id_2")]
tree_data[, deltaYear := year_4 - year_2]

# Melt to long format
tree_data = melt(tree_data, measure = patterns("^dbh", "^year"), value.name = c("dbh", "year"))

## Subset to the species of interest
# Get species code
code = speciesInfo[stri_detect(Nombre, regex = species), Code]
if (any(speciesInfo[Code == code, Code2] != code | speciesInfo[Code == code, Code3] != code))
	warning("Not sure the good species was selected!")

# Subset
tree_data = tree_data[species_code == code]
mostCommon_deltaYear = tree_data[, .N, by = deltaYear][, deltaYear[N == max(N)]]
tree_data = tree_data[deltaYear == mostCommon_deltaYear]

tree_data[, c("species_code", "variable", "deltaYear") := NULL]

# Remove measurements posterior to 2018 (no climatic data)
tree_data = tree_data[year <= 2018]

if (all.equal(tree_data[, plot_id_4], paste0(tree_data[, plot_id_2], "A1")))
{
	print("plot_id_4 and plot_id_2 are equivalent. I just keep plot_id_2 and renamed the column plot_id")
	tree_data[, plot_id_4 := NULL]
	setnames(x = tree_data, old = "plot_id_2", new = "plot_id")
} else {
	stop("I assumed that plot_id_4 and plot_id_2 are equivalent in the rest of this program!")
}

setkey(tree_data, plot_id, tree_id, year)
setcolorder(tree_data)

if (any(tree_data[, .N, by = .(plot_id, tree_id)][, N] < 2))
	stop("Some trees have only one measurement")

## Compute growth
growth_dt = copy(tree_data)
computeDiametralGrowth(growth_dt, byCols = c("plot_id", "tree_id"))
growth_dt = na.omit(growth_dt)

extreme = quantile(growth_dt[growth >= 0, growth], seq(0, 1, 0.001))["99.9%"]
remarquable_trees = unique(growth_dt[(growth < 0) | (growth > extreme), .(plot_id, tree_id)])

#### Get the environment time series (climate, soil pH, basal area) for this dataset
## Function to get all the years within a time range
getYears = function(time_range)
{
	years = stri_split(str = time_range, regex = "-", simplify = TRUE)
	year_1 = as.integer(years[1])
	year_2 = as.integer(years[2])

	return(as.character((year_1:year_2)))
}

## Get time-plot_id coordinates
spanishCoords = unique(tree_data[, .(plot_id, lon, lat)])
index_min = unique(tree_data[tree_data[, .(m = .I[year == min(year)]), by = .(plot_id, lon, lat)]$m,
	.(plot_id, lon, lat, year)])
index_max = unique(tree_data[tree_data[, .(M = .I[year == max(year)]), by = .(plot_id, lon, lat)]$M,
	.(plot_id, lon, lat, year)])

spanishCoords[index_min, on = "plot_id", min_year := i.year]
spanishCoords[index_max, on = "plot_id", max_year := i.year]
spanishCoords[, unique_year_id := paste(min_year, max_year, sep = "-")]

## Common variables
selectedVariables = c("pr", "tas")
clim_folder = "/bigdata/Predictors/climateChelsa/"
soil_folder = "/bigdata/Predictors/Soil\ esdacph\ Soil\ pH\ in\ Europe/"

## Compute length vector (sum of the number of years required for each plot)
length_clim = spanishCoords[, sum(max_year - min_year + 1)]
combination_years = unique(spanishCoords[, unique_year_id])

clim_dt = data.table()

for (j in selectedVariables)
	set(clim_dt, i = NULL, j = j, value = rep(0, length_clim))

clim_dt[, plot_id := "0"]
clim_dt[, year := 0]

count_clim_var = 0

## Stack raster climate
treeYears = tree_data[, min(year)]:tree_data[, max(year)]
for (clim_var in selectedVariables)
{
	clim_rs = rast(x = paste0(clim_folder, clim_var, "/", treeYears, ".tif"))
	names(clim_rs) = as.character(treeYears)
	count = 0
	start_ind = 1
	
	# Extraction
	for (time_range_id in combination_years)
	{
		count = count + 1
		selectedCoords = vect(spanishCoords[unique_year_id == time_range_id, .(lon, lat)], crs = "EPSG:23030",
			geom = c("lon", "lat")) # Order MUST be longitude, latitude
		selectedCoords = project(selectedCoords, "EPSG:4326")
		selectedPoint_id = spanishCoords[unique_year_id == time_range_id, plot_id]
		selectedLayers = getYears(time_range_id)
		
		values = setDT(extract(x = terra::subset(x = clim_rs, subset = selectedLayers), y = selectedCoords))
		values[, plot_id := selectedPoint_id]
		values = data.table::melt(data = values, id.vars = c("plot_id"), measure.vars = selectedLayers,
			variable.name = "year", value.name = "climate", variable.factor = FALSE)
		values[, year := as.integer(year)]
		end_ind = start_ind + values[, .N] - 1
		clim_dt[start_ind:end_ind, c("plot_id", "year") := values[, .(plot_id, year)]]
		clim_dt[start_ind:end_ind, c(clim_var) := values[, climate]]

		print(paste0(round(count*100/length(combination_years)), "% done for variable ", clim_var))
		start_ind = end_ind + 1
	}
	count_clim_var = count_clim_var + 1
	print(paste0("--->", round(count_clim_var*100/length(selectedVariables)), "% total done"))
}

setkey(clim_dt, plot_id, year)
setcolorder(clim_dt)

## Soil data. According to https://esdac.jrc.ec.europa.eu/content/soil-ph-europe, the projection is ETRS89 Lambert Azimuthal Equal Area 
soil = rast(paste0(soil_folder, "ph_cacl2/w001001.adf"))
soil_shp = vect(paste0(soil_folder, "country_laea.shp"))
crs(soil) = crs(soil_shp)

soil = project(x = soil, y = "EPSG:4326")

selectedCoords = vect(spanishCoords[, .(plot_id, lon, lat)], crs = "EPSG:23030",
	geom = c("lon", "lat")) # Order MUST be longitude, latitude
selectedCoords = project(selectedCoords, "EPSG:4326")

soil_values = extract(x = soil, y = selectedCoords, xy = TRUE, cell = TRUE)
sum(is.na(soil_values))
setDT(soil_values)

# Few values are NA (the offset of some forest plots fall into the ocean, river, lakes...)
problematicPlots = soil_values[is.na(ph_cacl2)]

if (problematicPlots[, .N] != 0)
{
	for (i in 1:problematicPlots[, .N])
		problematicPlots[i, ph_cacl2 := mean(soil[adjacent(x = soil, cells = problematicPlots[i, cell])]$ph_cacl2, na.rm = TRUE)]

	sum(is.na(problematicPlots))
	soil_values[problematicPlots[, ID], ph_cacl2 := problematicPlots[, ph_cacl2]]
}

soil_values[, plot_id := selectedCoords$plot_id]
soil_values[, c("ID", "cell", "x", "y") := NULL]

setnames(soil_values, old = "ph_cacl2", new = "ph")
setkey(soil_values, plot_id)
setcolorder(soil_values)

## Merge climate with soil data
clim_dt[, plot_id := as.integer(plot_id)]
env = soil_values[clim_dt, on = "plot_id"]

setkey(env, plot_id, year)

setcolorder(env)

#### Run prediction based on the French data
## Comments
# The aim is to validate the fitted model. For this, I run predictions on the French data. The starting dbh has been estimated previously,
#	and the ending dbh is measured. The climate data are ready so that we can extract the time series of the environmental conditions, including
#	soil and stand basal area.
# From https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html
#	if you plan on making many calls to $generate_quantities() then the most efficient option is to pass the paths of the CmdStan CSV output
#	files (this avoids CmdStanR having to rewrite the draws contained in the fitted model object to CSV each time). If you no longer have the CSV
#	files you can use draws_to_csv() once to write them and then pass the resulting file paths to $generate_quantities() as many times as needed.

## Load stand-alone generate_quantities() stan model
gq_model_ssm = cmdstan_model("generate_posterior.stan")

## Load results
info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

n_chains = ssm$num_chains()

## Prepare stan data
# Load
stanData_ssm = readRDS(paste0(tree_path, run, "_stanData.rds"))

# Add dimensions
stanData_ssm$n_climate_new = mostCommon_deltaYear + 1 # Number of latent dbh per individual
stanData_ssm$n_years = mostCommon_deltaYear # Number of growing years
stanData_ssm$n_growth = stanData_ssm$n_children
stanData_ssm$n_indiv_new = unique(tree_data[, .(plot_id, tree_id)])[, .N]

# Add diameters
parents_index = tree_data[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
stanData_ssm$dbh0 = tree_data[parents_index, dbh]

# Format and organise climate data
envStan_dt = data.table()

for (j in c("pr", "tas", "ph", "standBasalArea")) # The order matter, I assume it is later precipitations, temperatures, ph, basalAreas
	for (i in seq_len(stanData_ssm$n_climate_new))
		set(envStan_dt, i = NULL, j = paste0(j, i), value = rep(0, stanData_ssm$n_indiv_new))

count = 0

if (!all.equal(seq(1, tree_data[, .N], by = 2), parents_index))
	stop("I assumed that each tree is measured only twice!")

for (i in seq(1, tree_data[, .N], 2)) # Works only if tree_data sort by year for each tree within plot! Check keys = "plot_id" "tree_id" "year"
{
	# --- Common variables
	time_span = tree_data[i, year]:tree_data[i + 1, year]
	currentPlot = tree_data[i, plot_id]
	count = count + 1

	# --- Add the new environment (x_r)
	# ------ Approach: ssm
	precipitations = (env[.(currentPlot, time_span), pr] - stanData_ssm$pr_mu)/stanData_ssm$pr_sd
	temperatures = (env[.(currentPlot, time_span), tas] - stanData_ssm$tas_mu)/stanData_ssm$tas_sd
	ph = (env[.(currentPlot, time_span), ph] - stanData_ssm$ph_mu)/stanData_ssm$ph_sd
	basalAreas = rep(0, mostCommon_deltaYear + 1)

	envStan_dt[count, names(envStan_dt) := as.list(c(precipitations, temperatures, ph, basalAreas))]

	if (count %% 1000 == 0)
		print(paste0(round(i*100/tree_data[, .N], 2), "% done"))
}

stanData_ssm$x_r = envStan_dt

# cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "pr[:digit:]")]
# precipitations = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

# cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "tas[:digit:]")]
# temperatures = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

# cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "ph[:digit:]")]
# ph = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

# cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "standBasalArea[:digit:]")]
# basalAreas = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

# stanData_classic$x_r_avg = data.table(pr = precipitations, tas = temperatures, ph = ph, standBasalArea = basalAreas)

## Run simulations
# SSM aproach
generate_quantities_ssm = gq_model_ssm$generate_quantities(ssm$draws(), data = stanData_ssm, parallel_chains = n_chains)
simulated_dbh_ssm = stanData_ssm$sd_dbh*generate_quantities_ssm$draws("current_dbh")
increment_deltaYears_ssm = simulated_dbh_ssm[, , paste0("current_dbh[", 1:(stanData_ssm$n_indiv_new), ",", stanData_ssm$n_climate_new, "]")] - 
	simulated_dbh_ssm[, , paste0("current_dbh[", 1:(stanData_ssm$n_indiv_new), ",1]")]

residuals_ssm = increment_deltaYears_ssm
measuredIncrements = unique(growth_dt[, .(plot_id, tree_id, increment)])
for (i in 1:measuredIncrements[, .N]) # Less than a minute run for 46000 individuals
{
	residuals_ssm[, , i] = increment_deltaYears_ssm[, , i] - measuredIncrements[i, increment]
	if (i %% 100 == 0)
		print(paste0(round(i*100/measuredIncrements[, .N], 2), "% done"))
	# residuals_classic[, , i] = increment_5yrs_classic[, , i] - measuredIncrements[i, increment]
}

final_dbh = apply(X = simulated_dbh_ssm[, , paste0("current_dbh[", 1:(stanData_ssm$n_indiv_new), ",", stanData_ssm$n_climate_new, "]")],
	MARGIN = 3, FUN = mean)

growth_dt[, final_dbh_obs := dbh + increment]
growth_dt[, final_dbh_est := final_dbh]

# # Classic approach
# # generate_quantities_classic = gq_model_classic$generate_quantities(classic$draws(), data = stanData_classic, parallel_chains = n_chains)
# # simulated_dbh_classic = stanData_classic$sd_dbh*generate_quantities_classic$draws("current_dbh")
# # increment_5yrs_classic = simulated_dbh_classic[, , paste0("current_dbh[", 1:(stanData_classic$n_indiv_new), ",", stanData_ssm$n_years, "]")] - 
# # 	simulated_dbh_classic[, , paste0("current_dbh[", 1:(stanData_classic$n_indiv_new), ",1]")]

# ## Compute residuals
# measuredIncrements = unique(tree_data[, .(plot_id, tree_id, radial_increment_5_yrs_mm)])
# if (measuredIncrements[, .N] != dim(increment_5yrs_classic)[3])
# 	stop("Dimensions mismatch between measured vs simulated increments")

# residuals_ssm = increment_5yrs_ssm
# # residuals_classic = increment_5yrs_classic

# residuals_mean_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = mean)
# # residuals_mean_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = mean)

# residuals_median_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = median)
# # residuals_median_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = median)

# residuals_annualIncrement = apply(X = simulatedAnnualIncrement, MARGIN = 3, FUN = mean) - 2*measuredIncrements[, radial_increment_5_yrs_mm]/5

# for (i in 1:measuredIncrements[, .N]) # Less than a minute run for 46000 individuals
# {
# 	residuals_ssm[, , i] = increment_5yrs_ssm[, , i] - measuredIncrements[i, radial_increment_5_yrs_mm]
# 	# residuals_classic[, , i] = increment_5yrs_classic[, , i] - measuredIncrements[i, radial_increment_5_yrs_mm]
# }

##! This is the stuff I used to plot
# This reproduces the figure I sent to Lisa Friday the 27th January, at 17:12, for Fagus sylvatica. It is also in Fagus' folder on $local
# lazyComparePosterior(list(residuals_ssm[, , 3], residuals_classic[, , 3]), legend = c("ssm", "classic"))

# lazyComparePosterior(list(residuals_ssm, residuals_classic), legend = c("ssm", "classic"))
lazyPosterior(residuals_ssm, legend = "ssm")

colours = MetBrewer::met.brewer("Austria", 2)
colours_str = grDevices::colorRampPalette(colours)(2)
colours_str_pol = paste0(colours_str, "66")

ind = which(abs(residuals_ssm) < 200)
aa = density(residuals_ssm[ind], n = 512)

mean(residuals_ssm[ind])

plot(0, type = "n", xlim = c(min(aa$x), max(aa$x)), ylim = c(0, max(aa$y)), ylab = "frequence", xlab = "", main = "Posteriors")
lines(aa, col = colours_str[1], lwd = 4, lty = 1)
polygon(aa, col = colours_str_pol[1])
legend(x = "topright", legend = c("ssm mean", "classic mean"), fill = colours_str_pol, box.lwd = 0)



##! -- END This is the stuff I used to plot


##! Plot resiudals vs observed increments
residuals_ssm_mean = apply(X = residuals_ssm, MARGIN = 3, FUN = mean)
lm_residuals = lm(residuals_ssm_mean ~ growth_dt[, dbh])
lm_title = paste(round(lm_residuals$coefficients["(Intercept)"], 2), "+", round(lm_residuals$coefficients["growth_dt[, dbh]"], 2), "dbh")

plot(growth_dt[, dbh], residuals_ssm_mean, pch = 19, xlab = "Observed dbh", ylab = "Residuals predicted - observed dbh 'children'",
	main = lm_title, col = "#88779933")
abline(lm_residuals, lwd = 4, col = "#CD212A")
growth_dt[, ind := .I]
ind_weird = growth_dt[remarquable_trees, ind]
points(growth_dt[ind_weird, dbh], residuals_ssm_mean[ind_weird], pch = 19,
	col = "#CD212A")


residuals_ssm_mean_noWeirdos = residuals_ssm_mean[-ind_weird]
lm_residuals_noWeirdos = lm(residuals_ssm_mean_noWeirdos ~ growth_dt[-ind_weird, dbh])
lm_title = paste(round(lm_residuals$coefficients["(Intercept)"], 2), "+", round(lm_residuals$coefficients["growth_dt[, dbh]"], 2), "dbh")

plot(growth_dt[-ind_weird, dbh], residuals_ssm_mean_noWeirdos, pch = 19, xlab = "Observed dbh", ylab = "Residuals predicted - observed dbh 'children'",
	main = lm_title, col = "#88779933")
abline(lm_residuals, lwd = 4, col = "#CD212A")


##! Mean and median stuff

residuals_mean_ssm = residuals_mean_ssm - measuredIncrements[, radial_increment_5_yrs_mm]
residuals_mean_classic = residuals_mean_classic - measuredIncrements[, radial_increment_5_yrs_mm]
residuals_median_ssm = residuals_median_ssm - measuredIncrements[, radial_increment_5_yrs_mm]
residuals_median_classic = residuals_median_classic - measuredIncrements[, radial_increment_5_yrs_mm]

##! Plot residuals vs observed
plot(measuredIncrements[, radial_increment_5_yrs_mm], residuals_mean_ssm, xlab = "Observed increment",
	ylab = "Residuals on mean", pch = 19)
abline(lm(residuals_mean_ssm ~ measuredIncrements[, radial_increment_5_yrs_mm]), lwd = 2, col = "#CD212A")
abline(h = 0, lwd = 2, col = "#887799")
plot(measuredIncrements[, radial_increment_5_yrs_mm], residuals_median_ssm, xlab = "Observed increment",
	ylab = "Residuals on median", pch = 19)

##! Plot predictions vs observed
prediction_mean_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = mean)
prediction_mean_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = mean)

prediction_median_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = median)
prediction_median_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = median)

plot(measuredIncrements[, radial_increment_5_yrs_mm], prediction_mean_ssm, xlab = "Observed increment",
	ylab = "Mean predicted increment", pch = 19)
abline(a = 0, b = 1, lwd = 2, col = "#CD212A")
plot(measuredIncrements[, radial_increment_5_yrs_mm], prediction_median_ssm, xlab = "Observed increment",
	ylab = "Median predicted increment", pch = 19)
abline(a = 0, b = 1, lwd = 2, col = "#CD212A")

##! Plot on mean
colours = MetBrewer::met.brewer("Austria", 2)
colours_str = grDevices::colorRampPalette(colours)(2)
colours_str_pol = paste0(colours_str, "66")

density_ls = vector(mode = "list", length = 2)
density_ls[[1]] = density(residuals_mean_ssm, n = 512)
density_ls[[2]] = density(residuals_mean_classic, n = 512)

min_x = min(density_ls[[1]]$x, density_ls[[2]]$x)
max_x = max(density_ls[[1]]$x, density_ls[[2]]$x)
max_y = max(density_ls[[1]]$y, density_ls[[2]]$y)

plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", xlab = "", main = "Posteriors on mean")
lines(density_ls[[1]], col = colours_str[1], lwd = 4, lty = 1)
lines(density_ls[[2]], col = colours_str[2], lwd = 4, lty = 2)
polygon(density_ls[[1]], col = colours_str_pol[1])
polygon(density_ls[[2]], col = colours_str_pol[2])
legend(x = "topright", legend = c("ssm mean", "classic mean"), fill = colours_str_pol, box.lwd = 0)

##! Plot on median
density_ls = vector(mode = "list", length = 2)
density_ls[[1]] = density(residuals_median_ssm, n = 512)
density_ls[[2]] = density(residuals_median_classic, n = 512)

min_x = min(density_ls[[1]]$x)
max_x = max(density_ls[[1]]$x)
max_y = max(density_ls[[1]]$y)

plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", xlab = "", main = "Posteriors on median")
lines(density_ls[[1]], col = colours_str[1], lwd = 4, lty = 1)
polygon(density_ls[[1]], col = colours_str_pol[1])
legend(x = "topright", legend = c("ssm median"), fill = colours_str_pol, box.lwd = 0)

##! Plot on simulated observations (mean)
density_ls = vector(mode = "list", length = 1)
density_ls[[1]] = density(residuals_annualIncrement, n = 512)

min_x = min(density_ls[[1]]$x, density_ls[[2]]$x)
max_x = max(density_ls[[1]]$x, density_ls[[2]]$x)
max_y = max(density_ls[[1]]$y, density_ls[[2]]$y)

plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", xlab = "", main = "Posteriors on median")
lines(density_ls[[1]], col = colours_str[1], lwd = 4, lty = 1)
lines(density_ls[[2]], col = colours_str[2], lwd = 4, lty = 2)
polygon(density_ls[[1]], col = colours_str_pol[1])
polygon(density_ls[[2]], col = colours_str_pol[2])
legend(x = "topright", legend = c("ssm median", "classic median"), fill = colours_str_pol, box.lwd = 0)

##! Others
mean(residuals_ssm)
mean(residuals_classic)

sd(residuals_ssm)
sd(residuals_classic)

##! Others
obsPositiveAnnualAvgGrowth = stanData_ssm$avg_yearly_growth_obs[stanData_ssm$avg_yearly_growth_obs >= 0]
hist(obsPositiveAnnualAvgGrowth, xlab = "Observed positive annual averaged growth", breaks = seq(0, 35, 1))
hist(measuredIncrements[, 2*radial_increment_5_yrs_mm/5], xlab = "Observed averaged increment over 5 years", breaks = seq(0, 35, 1))

mean(obsPositiveAnnualAvgGrowth)
mean(measuredIncrements[, 2*radial_increment_5_yrs_mm/5])

sd(obsPositiveAnnualAvgGrowth)
sd(measuredIncrements[, 2*radial_increment_5_yrs_mm/5])

median(obsPositiveAnnualAvgGrowth)
median(measuredIncrements[, 2*radial_increment_5_yrs_mm/5])
