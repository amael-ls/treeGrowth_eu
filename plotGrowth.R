
#### Aim of prog: Plot different visualisations of growth versus climate and competition

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(terra)

#### Tool functions
## Growth function
growth_fct = function(x, precipitation, temperature, basalArea, params, scaled_dbh = FALSE, scaled_clim = FALSE, ...)
{
	providedArgs = list(...)
	providedArgs_names = names(providedArgs)
	scaling = 1

	if (!scaled_dbh)
	{
		if (!("sd_dbh" %in% providedArgs_names))
			stop("You need to provide sd_dbh in order to standardise dbh")
		
		sd_dbh = providedArgs[["sd_dbh"]]
		scaling = sd_dbh
		x = x/sd_dbh
	}

	if (!scaled_clim)
	{
		if (!(c("pr_mu", "pr_sd", "tas_mu", "tas_sd") %in% providedArgs_names))
			stop("You need to provide pr_mu, pr_sd, tas_mu, and tas_sd in order to standardise dbh")
		
		temperature = (temperature - tas_mu)/tas_sd
		precipitation = (precipitation - pr_mu)/pr_sd
	}

	potentialGrowth = params["potentialGrowth"]
	dbh_slope = params["dbh_slope"]

	pr_slope = params["pr_slope"]
	pr_slope2 = params["pr_slope2"]
	tas_slope = params["tas_slope"]
	tas_slope2 = params["tas_slope2"]

	competition_slope = params["competition_slope"]

	return(scaling*exp(potentialGrowth + dbh_slope*x + pr_slope*precipitation + pr_slope2*precipitation^2 +
		tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*basalRea))
}

#### Load data
## Common variables
species = "Tilia_platyphyllos"
year = 2015

## Paths
tree_path = paste0("./", species, "/")
climate_path = "/home/amael/project_ssm/climateData/Chelsa/yearlyAverage/"
shapefile_path = "/home/amael/shapefiles/deutschland/"
# shapefile_path = "/Users/mistral/Nextcloud/shapefiles/Germany_shapefile/"

## Tree data and parameters
results = readRDS(paste0(tree_path, "growth-2022-03-10_17h33.rds"))

## Climate
temperature_file = paste0(climate_path, "tas/", year, ".tif")
precipitation_file = paste0(climate_path, "pr/", year, ".tif")

climate = rast(c(temperature_file, precipitation_file))

## Landscape
germany = vect(paste0(shapefile_path, "germany.shp"))
germany = simplifyGeom(x = germany, tolerance = 100)
germany = project(x = germany, y = climate)

## Crop climate to Germany only
climate = crop(x = climate, y = germany, mask = TRUE)

## Scaling
dbh_mu_sd = readRDS(paste0(tree_path, "dbh_normalisation.rds"))
climate_mu_sd = readRDS(paste0(tree_path, "climate_normalisation.rds"))
