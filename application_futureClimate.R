
#### Aim of prog: Compute the probability G_2050 < median 2000-2010 for Fagus sylvatica
## Explanations
#	Averaged climatic data are used to compute the median of growth for each pixel. Then, the aim is to compute the probability
#		that growth in 2050 is below this median. If it is more than 50 %, then it means that there is a shift/contraction of the
#		distribution to the left, i.e., growth is likely to decrease by 2050
#
#	We only use Fagus sylvatica because it responds to climatic variables and because it is an important tree in west Europe
#
# 	The climatic variables and units for the future climate are:
#		bio1	= Annual Mean Temperature [°C*10]
#		bio12	= Annual Precipitation [mm/year]
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
## Source functions
source("toolFunctions.R")

## Compute meanlog posterior
gq_proba = function(climate, stanData, draws, gq_model, n_chains, N_draws, dbh0)
{
	stanData$N_draws = N_draws
	stanData$dbh0 = dbh0
	stanData$n_growth = stanData$n_obs - stanData$n_indiv
	stanData$env_current = c(climate[c("pr_current",  "tas_current")], 0, 0)
	stanData$env_future = c(climate[c("pr_future",  "tas_future")], 0, 0)

	if (any(is.na(stanData$env_current)) || any(is.na(stanData$env_future)))
		return (NA)
	
	generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

	return (mean(generate_quantities$draws("probaGrowth_inf_median")))
}

# aa = generate_quantities$draws("median_current")
# bb = generate_quantities$draws("median_future")

# stanData$env_current[1]*scaling["pr", sd] + scaling["pr", mu]
# unname(stanData$env_future[1])*scaling["pr", sd] + scaling["pr", mu]

# stanData$env_current[2]*scaling["tas", sd] + scaling["tas", mu]
# unname(stanData$env_future[2])*scaling["tas", sd] + scaling["tas", mu]


# qq = terra::app(x = test, fun = gq_proba, cores = 1, stanData = stanData, draws = draws,
# 	gq_model = gq_model, n_chains = n_chains, N_draws = 100, dbh0 = 0)

#### Prepare climatic data
## Load Fagus sylvatica distribution
fagus_distrib = vect("/home/amael/shapefiles/trees/europe/chorological_maps/Fagus sylvatica/shapefiles/Fagus_sylvatica_plg.shp")

## Load climate, average, and crop (i.e., mask from Fagus sylvatica distribution)
pr_2000_2010 = avgAnnualClimate(2000, 2010, "pr", crop = TRUE, area = fagus_distrib)
tas_2000_2010 = avgAnnualClimate(2000, 2010, "tas", crop = TRUE, area = fagus_distrib)

## Set the values of the distribution that falls into the sea/ocean
# Load europe shapefile
europe = vect("/home/amael/shapefiles/europe/continent/europe.shp")
europe = project(europe, pr_2000_2010)

# Masking and set NA
pr_2000_2010 = mask(x = pr_2000_2010, mask = europe)
tas_2000_2010 = mask(x = tas_2000_2010, mask = europe)

## Load and prepate climate 2050
# Folder
climate_folder = "~/scratch/Chelsa/"
if (!dir.exists(climate_folder))
	stop(paste0("The climate folder <", climate_folder, "> does not exist! Download the CMIP6 data first"))

variables = c("bio1", "bio12")
clim_model = "mpi-esm1-2-hr"

# Loading
pr_2050 = rast(paste0(climate_folder, variables[2], "/", clim_model, "/", variables[2], "_", clim_model, ".tif"))
tas_2050 = rast(paste0(climate_folder, variables[1], "/", clim_model, "/", variables[1], "_", clim_model, ".tif"))

# Projecting, croping, and masking
pr_2050 = project(x = pr_2050, y = pr_2000_2010)
tas_2050 = project(x = tas_2050, y = tas_2000_2010)

pr_2050 = crop(x = pr_2050, y = pr_2050)
tas_2050 = crop(x = tas_2050, y = tas_2050)

pr_2050 = mask(x = pr_2050, mask = fagus_distrib)
tas_2050 = mask(x = tas_2050, mask = fagus_distrib)

climate_rs = c(pr_2000_2010, tas_2000_2010, pr_2050, tas_2050)
names(climate_rs) = c("pr_current", "tas_current", "pr_future", "tas_future")

#### Compute proba future growth < current median
## Load results to get posterior parameters
species = "Fagus sylvatica"
path = paste0("./", species, "/")
run = 1
lastRun = getLastRun(path = path, run = run)

results = readRDS(paste0(path, lastRun[["file"]]))

## Prepare data for generated quantities model to run
stanData = readRDS(paste0(path, run, "_stanData.rds"))
n_chains = results$num_chains()
draws = results$draws()

scaling = readRDS(paste0(path, run, "_climate_normalisation.rds"))
setkey(scaling, variable)

## Compile model for posterior median, gq = generated quantities
gq_model = cmdstan_model("./generate_probaLower.stan")

coords = as.polygons(ext(10, 10.04, 50, 50.04), crs = crs(climate_rs))
# coords = ext(10, 11, 50, 51)
test = terra::extract(x = climate_rs, y = coords, xy = TRUE)
setDT(test)
setcolorder(test, c("x", "y"))
test[, ID := NULL]
test = rast(test, type = "xyz")

test[["pr_current"]] = (test[["pr_current"]] - scaling["pr", mu])/scaling["pr", sd]
test[["tas_current"]] = (test[["tas_current"]] - scaling["tas", mu])/scaling["tas", sd]

test[["pr_future"]] = (test[["pr_future"]] - scaling["pr", mu])/scaling["pr", sd]
test[["tas_future"]] = (test[["tas_future"]] - scaling["tas", mu])/scaling["tas", sd]

## Compute probabilities
start = Sys.time()
toto = terra::app(x = test, fun = gq_proba, cores = 1, stanData = stanData, draws = draws,
	gq_model = gq_model, n_chains = n_chains, N_draws = 100, dbh0 = 0)
end = Sys.time()
