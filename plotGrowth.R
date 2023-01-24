
#### Aim of prog: Plot different visualisations of growth versus predictors (diameters, environmental factors)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
source("toolFunctions.R")

#### Load data
## Common variables
args = commandArgs(trailingOnly = TRUE)
# args = c("Fagus sylvatica", "1", "pr")
if (length(args) != 3)
	stop("Supply the species_id, run_id, and max_indiv as command line arguments!", call. = FALSE)

species = as.character(args[1])
run = as.integer(args[2])
variable = as.character(args[3])

if (!(variable %in% c("pr", "tas", "ph")))
	stop("The variable must be either pr, tas, or ph")

info = plotGrowth(species, run, variable)
print(info)

########! START CRASH TEST ZONE -----------------------------------------------------------------------------------------------------------------
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
########! END CRASH TEST ZONE -------------------------------------------------------------------------------------------------------------------