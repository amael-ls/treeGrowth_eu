
#### Aim of prog: Plot different visualisations of growth versus predictors (diameters, environmental factors)
## Comments
#	The first part is dedicated to the growth response to environment and to time series predictions
#	The second part is dedicated to 'validation'
#		For this part, I use radial increment based on wood cores from Martínez-Sancho 2020: https://doi.org/10.1038/s41597-019-0340-y
#			Dendroecological Collection, tree-ring and wood density data from seven tree species across Europe.
#		These data have not been used to parametrise the growth model because they do not match the error structure used in our model
#		these data are only available for the 7 following species:
#			- Betula pendula
#			- Fagus sylvatica
#			- Picea abies
#			- Pinus pinaster
#			- Pinus sylvestris
#			- Populus nigra
#			- Quercus petraea
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

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)
library(terra)

#### Tool functions
source("toolFunctions.R")

######## Part I: growth response to environment and to time series predictions
#### Get data
## Common variables
# args = c("Betula pendula", "1", "tas", "pr", "ph", "ba")
# args = c("Fagus sylvatica", "1", "tas", "pr", "ph", "ba")
# args = c("Picea abies", "1", "tas", "pr", "ph", "ba")
# args = c("Pinus pinaster", "1", "tas", "pr", "ph", "ba")
# args = c("Pinus sylvestris", "1", "tas", "pr", "ph", "ba")
# args = c("Quercus petraea", "1", "tas", "pr", "ph", "ba")
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
	stop("Supply the species_id, run_id, and at least one variable among pr, tas, and ph!", call. = FALSE)

species = as.character(args[1])
run = as.integer(args[2])
variable = as.character(args[3:length(args)])

if (!all(variable %in% c("pr", "tas", "ph", "ba")))
	stop("The variable must be either pr, tas, ph, or ba")

#### Plot residuals fit
# residuals_fit(species, run, filenamePattern = "_residuals_fit.pdf")

#### Plot growth vs variable
if ((file.exists("./speciesInformations.rds")) && (file.exists("./speciesInformations_runs.rds")))
{
	info = readRDS("./speciesInformations.rds")
	info_runs = readRDS("./speciesInformations_runs.rds")
	ls_info = list(info = info, range_subset = info_runs)
} else {
	ls_info = infoSpecies()
}

ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
info_runs = info_runs[ls_species] # speciesName_sci is the key of info_runs, so no need for %in% ls_species
info = info[ls_species] # speciesName_sci is the key of info, so no need for %in% ls_species

# All the variables range are within the species range, i.e., min(xyz_025) > min(xyz) among all species, and same with max
checkup = plotGrowth(species = species, run = run, ls_info = ls_info, variables = variable, extension = "tex", caption = FALSE,
	span = list(pr = info[, c(min(pr_025), max(pr_975))], tas = info[, c(min(tas_025), max(tas_975))],
		ph = info[, c(min(ph_025), max(ph_975))], ba = info[, c(min(ba_025), max(ba_975))]),
	pr_min = info[species, pr_025], pr_max = info[species, pr_975],
	tas_min = info[species, tas_025], tas_max = info[species, tas_975],
	ph_min = info[species, ph_025], ph_max = info[species, ph_975],
	ba_min = info[species, ba_025], ba_max = info[species, ba_975])
# print(checkup)

plotGrowth_dbh(species = species, run = run, ls_info = ls_info, caption = FALSE, extension = "tex",
	dbh_min = info_runs[, min(dbh_025_1)], dbh_max = info_runs[, max(dbh_975_1)])


params_names = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope")


params_values = printParams(ls_species = ls_species, params_names = params_names, run = 1)

ssm_05 = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("q5"))
ssm_med = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("med"))
ssm_95 = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("q95"))

classic_05 = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("q5"))
classic_med = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("med"))
classic_95 = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("q95"))


data.table:::print.data.table(x = ssm_05, digits = 2)
data.table:::print.data.table(x = ssm_med, digits = 2)
data.table:::print.data.table(x = ssm_95, digits = 2)

# dataEnv_ls = getEnvSeries(species, run)

# validationTreeRing(species, run)

# speciesName_sci           plot_id tree_id   year   dbh
# Fagus sylvatica   germany_12591_3      10   2000   373
# Fagus sylvatica   germany_12591_3      10   2012   437

########! START CRASH TEST ZONE 0 ---------------------------------------------------------------------------------------------------------------
######## Part II: validation

# ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Populus nigra", "Quercus petraea")

# count = 0
# op = par(mar = c(2.5, 8, 0.1, 0.4), mgp = c(1.5, 0.4, 0), oma = c(0, 0, 0.9, 0), tck = -0.01, las = 1)
# plot(x = NULL, y = NULL, xlim = c(-1, 1), ylim = c(0.5, length(ls_species)), axes = FALSE, xlab = "Correlation", ylab = "")
# axis(side = 1, at = seq(-1, 1, by = 0.5))
# axis(side = 2, at = seq(1, length(ls_species), length.out = length(ls_species)),
# 	labels = ls_species, tck = 0, lwd = "")
# for (species in ls_species)
# {
# 	count = count + 1
# 	correl_dt = readRDS(paste0("./", species, "/dt_correl_spearman_fr-de-se.rds"))
# 	points(x = c(correl_dt[1:3, corr]), y = rep(count - 0.1, 3), pch = 18, col = c("#000000", "#11AED9", "#EF8A47"))
# 	points(x = c(correl_dt[7:8, corr]), y = rep(count, 2), pch = 19, col = c("#11AED9", "#EF8A47"))
# 	points(x = c(correl_dt[4:6, corr]), y = rep(count + 0.1, 3), pch = 15, col = c("#000000", "#11AED9", "#EF8A47"))
# 	abline(h = count, lwd = 0.5, col = "#AABBCC")
# }
# abline(v = 0, lwd = 0.5, col = "#445567")
# caption1 = legend(x = "bottomright", legend = c("Data", "SSM", "Classic", "Data", "Precip", "Temp"), pch = c(19, 18, 15, 22, 22, 22),
# 	col = c("#000000", "#000000", "#000000", "#000000", "#11AED9", "#EF8A47"), pt.bg = c("#000000", "#11AED9", "#EF8A47"),
# 	bty = "n", ncol = 2)

########! END CRASH TEST ZONE 0 -----------------------------------------------------------------------------------------------------------------



########! START CRASH TEST ZONE 1 ---------------------------------------------------------------------------------------------------------------
######## Part II: validation
#### Load data
## All data (used later to keep exclusively french data that are not involved in the parametrisation of growth)
treeFolder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(treeFolder))
	stop(paste0("Folder\n\t", treeFolder, "\ndoes not exist"))

treeData = readRDS(paste0(treeFolder, "standardised_european_growth_data_reshaped.rds"))
setkey(treeData, plot_id, tree_id, year)
ls_species = sort(treeData[, unique(speciesName_sci)])

## French data
frenchFolder = "/bigdata/Inventories/FR IFN/processed data/"
if (!dir.exists(frenchFolder))
	stop(paste0("The folder <", frenchFolder, "> does not exist!"))

frenchData = readRDS(paste0(frenchFolder, "trees_forest.rds"))
stand_BA_dt = unique(readRDS(paste0(frenchFolder, "trees_forest_reshaped.rds"))[, .(pointInventory_id, year, standBasalArea)])

#### Subset data
## Keep only columns of interest and adapt classes or primary key columns to treeData classes
keptCols = c("year", "tree_id", "pointInventory_id", "speciesName_sci", "dbh_in_mm", "state_1stVisit", "radial_increment_5_yrs_mm",
	"xLambert93", "yLambert93")
frenchData = frenchData[, ..keptCols]

setnames(frenchData, old = "pointInventory_id", new = "plot_id")

set(frenchData, j = "plot_id", value = as.character(frenchData[["plot_id"]]))
set(frenchData, j = "tree_id", value = as.character(frenchData[["tree_id"]]))
set(frenchData, j = "year", value = as.numeric(frenchData[["year"]]))

frenchData[, plot_id := paste0("france_", plot_id)]
frenchData = frenchData[year <= 2018] # Maximum year available for climate data

setkey(frenchData, plot_id, tree_id, year)

## Adapt classes or primary key columns to treeData classes
setnames(stand_BA_dt, old = "pointInventory_id", new = "plot_id")
set(stand_BA_dt, j = "plot_id", value = as.character(stand_BA_dt[["plot_id"]]))
set(stand_BA_dt, j = "year", value = as.numeric(stand_BA_dt[["year"]]))

stand_BA_dt[, plot_id := paste0("france_", plot_id)]
stand_BA_dt = stand_BA_dt[year <= 2018] # Maximum year available for climate data

## Remove non-parametrised species
frenchData = frenchData[speciesName_sci %in% ls_species]

## Remove the NAs radial increment
frenchData = frenchData[!is.na(radial_increment_5_yrs_mm)]

## Keep only living trees
frenchData = frenchData[stri_detect(str = state_1stVisit, regex = "living")]

## Remove the remaining french data involved in the parametrisation of the growth model
frenchData = frenchData[!treeData]

## Check uniqueness individuals (I assume it so that I can greatly simplify the algorithm to estimate previous dbh)
uniquenessIndiv = all(frenchData[, .N, by = .(tree_id, plot_id)][, N] == 1)
if (!uniquenessIndiv)
	stop("I assumed that these trees have been measured only once!")

frenchData = frenchData[speciesName_sci == species]

#### Estimate the previous dbh
## Comments:
#	Because the last increment is not included, I need to estimate it. I assume the incomplete increment is 0.5 * incr/5
#		The 0.5 comes from the fact that I assume the tree is measured at the middle of the growing season, while incr/5
#		represents the averaged growth between t - 6 and t - 1, the year t being the current year for which the current
#		dbh is measured. At the end, incr + 0.5 * incr/5 = 11/10 incr.
#	Example: Imagine a tree is measured in June 2006. We can assume that the last tree ring is not complete yet, so the
#		increment over the last 5 growing years corresponds to the growth between 2000 and 2005:
#			1/ 2000 - 2001
#			2/ 2001 - 2002
#			3/ 2002 - 2003
#			4/ 2003 - 2004
#			5/ 2004 - 2005
#			6/ 2006: dbh is measured, but the tree ring of that year is not complete and is not accounted in increment_5yrs
#
#	Finally, because the increment is on the radius, I need to substract twice the increment to the diameter: dbh - 11/5 incr.
#

## Dbh at t - 6
# Compute
frenchData[, dbh_6years_ago := dbh_in_mm - 11*radial_increment_5_yrs_mm/5]
if (frenchData[, any(dbh_6years_ago < 0)])
{
	warning("There are some negative dbh!")
	print(frenchData[dbh_6years_ago < 0])
	frenchData = frenchData[dbh_6years_ago > 0]
}

frenchData[, years_minus6 := year - 6]

## Melting
frenchData = melt(frenchData, measure = patterns("^dbh", "^year"), value.name = c("dbh", "year"))
frenchData[, variable := NULL]

## Check that now each individual has two lines
if (!all(frenchData[, .N, by = .(tree_id, plot_id)][, N] == 2))
	stop("Check the computation of the previous dbh! Maybe 'melt' failed")

setkey(frenchData, plot_id, tree_id, year)

growth = (frenchData[seq(2, .N, by = 2), dbh] - frenchData[seq(1, .N, by = 2), dbh])/6
print(paste("Range of growth (in mm/year):", paste(round(range(growth), 2), collapse = " -> ")))

## Remove useless columns, reorder
frenchData[, c("state_1stVisit", "speciesName_sci") := NULL]
setcolorder(frenchData)

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
frenchCoords = unique(frenchData[, .(plot_id, xLambert93, yLambert93)])
index_min = unique(frenchData[frenchData[, .(m = .I[year == min(year)]), by = .(plot_id, xLambert93, yLambert93)]$m,
	.(plot_id, xLambert93, yLambert93, year)])
index_max = unique(frenchData[frenchData[, .(M = .I[year == max(year)]), by = .(plot_id, xLambert93, yLambert93)]$M,
	.(plot_id, xLambert93, yLambert93, year)])

frenchCoords[index_min, on = "plot_id", min_year := i.year]
frenchCoords[index_max, on = "plot_id", max_year := i.year]
frenchCoords[, unique_year_id := paste(min_year, max_year, sep = "-")]

## Common variables
selectedVariables = c("pr", "tas", "tasmin", "tasmax")
clim_folder = "/bigdata/Predictors/climateChelsa/"
soil_folder = "/bigdata/Predictors/Soil\ esdacph\ Soil\ pH\ in\ Europe/"

## Compute length vector (sum of the number of years required for each plot)
length_clim = frenchCoords[, sum(max_year - min_year + 1)]
combination_years = unique(frenchCoords[, unique_year_id])

clim_dt = data.table()

for (j in selectedVariables)
	set(clim_dt, i = NULL, j = j, value = rep(0, length_clim))

clim_dt[, plot_id := "0"]
clim_dt[, year := 0]

count_clim_var = 0

## Stack raster climate
treeYears = sort(unique(frenchData[, year]))
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
		selectedCoords = vect(frenchCoords[unique_year_id == time_range_id, .(xLambert93, yLambert93)], crs = "EPSG:2154",
			geom = c("xLambert93", "yLambert93")) # Order MUST be longitude, latitude
		selectedCoords = project(selectedCoords, "EPSG:4326")
		selectedPoint_id = frenchCoords[unique_year_id == time_range_id, plot_id]
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

selectedCoords = vect(frenchCoords[, .(plot_id, xLambert93, yLambert93)], crs = "EPSG:2154",
	geom = c("xLambert93", "yLambert93")) # Order MUST be longitude, latitude
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

## Merge climate with stand basal area and soil data
env = stand_BA_dt[clim_dt, on = c("plot_id", "year")]
env = env[soil_values, on = "plot_id"]

## Replace NAs in stand basal area by the unique non NA value of each plot
setkey(env, plot_id, year)

env[, isUniqueBA := sum(!is.na(standBasalArea)) == 1, by = plot_id]
if (any(env[, !isUniqueBA]))
	stop("I assumed that there is only one non NA basal area!")

env[, standBasalArea := env[plot_id][!is.na(standBasalArea), standBasalArea], by = plot_id]

env[, isUniqueBA := NULL]
env[, tasmin := NULL]
env[, tasmax := NULL]

#### Run prediction based on the French data
## Comments
# The aim is to validate the fitted model. For this, I run predictions on the French data. The starting dbh has been estimated previously,
#	and the ending dbh is measured. The climate data are ready so that we can extract the time series of the environmental conditions, including
#	soil and stand basal area.
# From https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html
#	if you plan on making many calls to $generate_quantities() then the most efficient option is to pass the paths of the CmdStan CSV output
#	files (this avoids CmdStanR having to rewrite the draws contained in the fitted model object to CSV each time). If you no longer have the CSV
#	files you can use draws_to_csv() once to write them and then pass the resulting file paths to $generate_quantities() as many times as needed.

## Common variables
tree_path = paste0("./", species, "/")
if (!dir.exists(tree_path))
	stop(paste0("Path not found for species <", species, ">."))

## Load stand-alone generate_quantities() stan models
gq_model_ssm = cmdstan_model("generate_posterior.stan")
gq_model_classic = cmdstan_model("generate_posterior_classic.stan")

## Load results
info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
classic = readRDS(paste0(tree_path, info_lastRun[["file"]]))

n_chains = classic$num_chains()

## Prepare stan data
# Load
stanData_ssm = readRDS(paste0(tree_path, run, "_stanData.rds"))
stanData_classic = readRDS(paste0(tree_path, run, "_stanData_classic.rds"))

# Add dimensions
stanData_ssm$n_climate_new = 7
stanData_ssm$n_years = 6
stanData_ssm$n_growth = stanData_ssm$n_children
stanData_ssm$n_indiv_new = unique(frenchData[, .(plot_id, tree_id)])[, .N]

stanData_classic$n_climate_new = 7
stanData_classic$n_years = 6
stanData_classic$n_indiv_new = stanData_ssm$n_indiv_new

# Add diameters
stanData_ssm$dbh0 = frenchData[seq(1, .N, by = 2), dbh]
stanData_classic$dbh0 = frenchData[seq(1, .N, by = 2), dbh]

# Format and organise climate data
envStan_dt = vector(mode = "list", length = 2)
names(envStan_dt) = c("ssm", "classic")

envStan_dt[["ssm"]] = data.table()
envStan_dt[["classic"]] = data.table()

for (modelType in names(envStan_dt))
	for (j in c("pr", "tas", "ph", "standBasalArea")) # The order matter, I assume it is later precipitations, temperatures, ph, basalAreas
		for (i in 1:stanData_classic$n_climate_new)
			set(envStan_dt[[modelType]], i = NULL, j = paste0(j, i), value = rep(0, stanData_ssm$n_indiv_new))

count = 0
for (i in seq(1, frenchData[, .N], 2))
{
	# --- Common variables
	time_span = frenchData[i, year]:frenchData[i + 1, year]
	currentPlot = frenchData[i, plot_id]
	count = count + 1

	# --- Add the new environment (x_r)
	# ------ Approach: ssm
	precipitations = (env[.(currentPlot, time_span), pr] - stanData_ssm$pr_mu)/stanData_ssm$pr_sd
	temperatures = (env[.(currentPlot, time_span), tas] - stanData_ssm$tas_mu)/stanData_ssm$tas_sd
	ph = (env[.(currentPlot, time_span), ph] - stanData_ssm$ph_mu)/stanData_ssm$ph_sd
	basalAreas = (env[.(currentPlot, time_span), standBasalArea] - stanData_ssm$ba_mu)/stanData_ssm$ba_sd

	envStan_dt[["ssm"]][count, names(envStan_dt[["ssm"]]) := as.list(c(precipitations, temperatures, ph, basalAreas))]

	# ------ Approach: classic
	precipitations = (env[.(currentPlot, time_span), pr] - stanData_classic$pr_mu)/stanData_classic$pr_sd
	temperatures = (env[.(currentPlot, time_span), tas] - stanData_classic$tas_mu)/stanData_classic$tas_sd
	ph = (env[.(currentPlot, time_span), ph] - stanData_classic$ph_mu)/stanData_classic$ph_sd
	basalAreas = (env[.(currentPlot, time_span), standBasalArea] - stanData_classic$ba_mu)/stanData_classic$ba_sd

	envStan_dt[["classic"]][count, names(envStan_dt[["classic"]]) := as.list(c(precipitations, temperatures, ph, basalAreas))]

	if (count %% 1000 == 0)
		print(paste0(round(i*100/frenchData[, .N], 2), "% done"))
}

stanData_ssm$x_r = envStan_dt[["ssm"]]
stanData_classic$x_r = envStan_dt[["classic"]]

cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "pr[:digit:]")]
precipitations = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "tas[:digit:]")]
temperatures = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "ph[:digit:]")]
ph = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

cols = names(envStan_dt[["classic"]])[stri_detect(names(envStan_dt[["classic"]]), regex = "standBasalArea[:digit:]")]
basalAreas = apply(X = envStan_dt[["classic"]][, ..cols], FUN = mean, MARGIN = 1)

stanData_classic$x_r_avg = data.table(pr = precipitations, tas = temperatures, ph = ph, standBasalArea = basalAreas)

## Run simulations
# SSM aproach
generate_quantities_ssm = gq_model_ssm$generate_quantities(ssm$draws(), data = stanData_ssm, parallel_chains = n_chains)
simulated_dbh_ssm = stanData_ssm$sd_dbh*generate_quantities_ssm$draws("current_dbh")
increment_5yrs_ssm = simulated_dbh_ssm[, , paste0("current_dbh[", 1:(stanData_ssm$n_indiv_new), ",", stanData_ssm$n_years, "]")] - 
	simulated_dbh_ssm[, , paste0("current_dbh[", 1:(stanData_ssm$n_indiv_new), ",1]")]

simulatedAnnualIncrement = stanData_ssm$sd_dbh*generate_quantities_ssm$draws("simulatedObservedGrowth")

# Classic approach
generate_quantities_classic = gq_model_classic$generate_quantities(classic$draws(), data = stanData_classic, parallel_chains = n_chains)
simulated_dbh_classic = stanData_classic$sd_dbh*generate_quantities_classic$draws("current_dbh")
increment_5yrs_classic = simulated_dbh_classic[, , paste0("current_dbh[", 1:(stanData_classic$n_indiv_new), ",", stanData_ssm$n_years, "]")] - 
	simulated_dbh_classic[, , paste0("current_dbh[", 1:(stanData_classic$n_indiv_new), ",1]")]

## Compute residuals
measuredIncrements = unique(frenchData[, .(plot_id, tree_id, radial_increment_5_yrs_mm)])
if (measuredIncrements[, .N] != dim(increment_5yrs_classic)[3])
	stop("Dimensions mismatch between measured vs simulated increments")

residuals_ssm = increment_5yrs_ssm
residuals_classic = increment_5yrs_classic

residuals_mean_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = mean)
residuals_mean_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = mean)

residuals_median_ssm = apply(X = increment_5yrs_ssm, MARGIN = 3, FUN = median)
residuals_median_classic = apply(X = increment_5yrs_classic, MARGIN = 3, FUN = median)

residuals_annualIncrement = apply(X = simulatedAnnualIncrement, MARGIN = 3, FUN = mean) - 2*measuredIncrements[, radial_increment_5_yrs_mm]/5

for (i in 1:measuredIncrements[, .N]) # Less than a minute run for 46000 individuals
{
	residuals_ssm[, , i] = increment_5yrs_ssm[, , i] - measuredIncrements[i, radial_increment_5_yrs_mm]
	residuals_classic[, , i] = increment_5yrs_classic[, , i] - measuredIncrements[i, radial_increment_5_yrs_mm]
}

##! This is the stuff I used to plot
# This reproduces the figure I sent to Lisa Friday the 27th January, at 17:12, for Fagus sylvatica. It is also in Fagus' folder on $local
# lazyComparePosterior(list(residuals_ssm[, , 3], residuals_classic[, , 3]), legend = c("ssm", "classic"))

lazyComparePosterior(list(residuals_ssm, residuals_classic), legend = c("ssm", "classic"))



ind = which(residuals_ssm < 50)
aa = density(residuals_ssm[ind], n = 512)

mean(residuals_ssm[ind])

plot(0, type = "n", xlim = c(min(aa$x), max(aa$x)), ylim = c(0, max(aa$y)), ylab = "frequence", xlab = "", main = "Posteriors")
lines(aa, col = colours_str[1], lwd = 4, lty = 1)
polygon(aa, col = colours_str_pol[1])
legend(x = "topright", legend = c("ssm mean", "classic mean"), fill = colours_str_pol, box.lwd = 0)



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

########! END CRASH TEST ZONE 1 -----------------------------------------------------------------------------------------------------------------

########! START CRASH TEST ZONE 2 ---------------------------------------------------------------------------------------------------------------
# #### Load data
# ## Common variables
# years = 2010:2018

# ## Paths
# climate_path = "/bigdata/Predictors/climateChelsa/"
# ph_path = "/bigdata/Predictors/Soil\ esdacph\ Soil\ pH\ in\ Europe/"
# shapefile_path = "/home/amael/shapefiles/Deutschland/"

# ## Climate (now, it would be more accurate to call this predictors, but by that time there were only climate!)
# temperature_files = paste0(climate_path, "tas/", years, ".tif")
# precipitation_files = paste0(climate_path, "pr/", years, ".tif")

# soil_shp = vect(paste0(ph_path, "country_laea.shp"))
# soil = rast(paste0(ph_path, "ph_cacl2/w001001.adf"))
# crs(soil) = crs(soil_shp)

# climate = rast(c(temperature_files, precipitation_files))
# soil = project(soil, climate)
# climate = c(climate, soil)
# names(climate) = c(paste0("tas_", years), paste0("pr_", years), "ph")

# ## Landscape
# germany = vect(paste0(shapefile_path, "germany.shp"))
# germany = simplifyGeom(x = germany, tolerance = 100)
# germany = project(x = germany, y = climate)

# ## Crop climate to Germany only
# climate = crop(x = climate, y = germany, mask = TRUE)
# climate_mean = c(mean(subset(climate, paste0("tas_", years))), mean(subset(climate, paste0("pr_", years))), mean(subset(climate, "ph")))
# names(climate_mean) = c("temperature", "precipitation", "ph")

# #### Compute average growth in a landscape/along a gradient
# ## Compute average growth in a landscape
# average_growth = setDT(as.data.frame(x = climate_mean, xy = TRUE))

# average_growth[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
# 	pr = precipitation, tas = temperature, ph = ph, basalArea = 25,
# 	params = estimated_params, standardised_dbh = FALSE, standardised_params = TRUE, standardised_variables = FALSE,
# 	sd_dbh = dbh_mu_sd[1, sd], scaling = rbindlist(list(climate_mu_sd, ph_mu_sd, ba_mu_sd)))$value,
# 	by = .(x, y)]

# average_growth[, integral := integral/(upper_dbh - lower_dbh)]

# #### Plot average growth
# ## In a landscape
# growth_rs = rast(x = average_growth, type = "xyz", crs = crs(climate_mean))

# png(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "average_growth_landscape.png"), width = 1080, height = 1080)
# plot(subset(growth_rs, "integral"), axes = FALSE, col = viridis(256),
# 	plg = list(cex = 3, title = "Growth(mm/yr)"))
# plot(germany, add = TRUE, col = NA, lwd = 2)
# plot(cities_rs, pch = 19, add = TRUE, cex = 4)
# dev.off()
########! END CRASH TEST ZONE 2 -----------------------------------------------------------------------------------------------------------------
