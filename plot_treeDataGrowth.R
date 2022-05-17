#### Aim of prog: Plot tree data on a map and few others informations (climate range, etc)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#? --------------------------------------------------------------------------------------------------------
######## PART I: Spatial plot of the data (plot location)
#? --------------------------------------------------------------------------------------------------------
#### Load data
## Folders
treeData_folder = "/home/amael/project_ssm/inventories/growth/"
shapefile_folder = "/home/amael/shapefiles/europe/"

## Tree data and shapefile
treeData = readRDS(paste0(treeData_folder, "standardised_european_growth_data_reshaped.rds"))
coords = vect(unique(treeData[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:4326")

europe = vect(paste0(shapefile_folder, "europe.shp"))

#### Plot
pdf("plots_location_growth.pdf", height = 10, width = 10)
plot(europe, col = "#C4AC7C44", border = "#9E391A", axes = FALSE)
plot(coords, pch = 20, cex = 0.025, col = "#354536", add = TRUE)
plot(europe, col = NA, border = "#9E391A", add = TRUE)
dev.off()



#? --------------------------------------------------------------------------------------------------------
######## PART II: Plot of the data in the climate and ph space, with gam models as a fit
#? --------------------------------------------------------------------------------------------------------
#### Tool function
## Create a colour gradient
gradColours = function(x, colours, n)
	return(colorRampPalette(colours) (n) [findInterval(x, seq(min(x), max(x), length.out = n))])

#### Common variables
## Folders
climate_folder = "/home/amael/project_ssm/inventories/growth/"
ph_folder = "/home/amael/project_ssm/inventories/growth/"

## Load climate and pH data
climate = readRDS(paste0(climate_folder, "europe_reshaped_climate.rds"))
ph = readRDS(paste0(ph_folder, "europe_reshaped_soil.rds"))

## Load species table
infoSpecies = readRDS("./speciesInformations.rds")

## Create explanatory data table for predictors
predictors_dt = data.table(predictors = c("pr", "tas", "tasmin", "tasmax"),
	explanation = c("annual precipitation", "mean annual temperature", "annual minimal temperature", "annual maximal temperature"))
predictors_dt[, complementary := ifelse(stri_detect(predictors, regex = "^tas"), "pr", "tas")]
setkey(predictors_dt, predictors)

## Load function
source("./extractClimate.R")

#### For loop for each species
for (i in 1:infoSpecies[, .N])
{
	species = infoSpecies[i, speciesName_sci]
	path = paste0("./", species, "/")
	if (!dir.exists(path))
		next;
	
	indices = readRDS(paste0(path, "full_indices.rds"))
	predictors_data = extractClimate(climate = climate, ph = ph, indices = indices)[["extracted_predictors"]]
	subset_trees = treeData[speciesName_sci == species, .(plot_id, tree_id, year, dbh)]

	# Compute growth. It is assumed that the data set is ordered by plot_id, tree_id, year!
	subset_trees[, prev_dbh := shift(dbh, 1, NA), by = .(plot_id, tree_id)]
	subset_trees[, prev_year := shift(year, 1, NA), by = .(plot_id, tree_id)]
	subset_trees[, annual_growth := (dbh - prev_dbh)/(year - prev_year)]

	# Cleaning data
	subset_trees[, c("prev_dbh", "prev_year") := NULL]
	subset_trees = subset_trees[!is.na(annual_growth)]

	predictors_data[, paste0(predictors_dt[, predictors], "_avg") := lapply(.SD, mean), .SDcols = predictors_dt[, predictors], by = plot_id]

	subset_trees = merge.data.table(subset_trees, predictors_data, by = c("plot_id", "year"))

	#### Plot
	## Set plot layout
	for (predictor in predictors_dt[, predictors])
	{
		pdf(paste0(path, "growth_vs_", predictor,".pdf"), height = 10, width = 10)
		layout(mat = matrix(c(2, 1, 0, 3, 0, 4), nrow = 2, ncol = 3),
			heights = c(1, 2), # Heights of the two rows
			widths = c(2, 1, 1)) # Widths of the three columns

		## Variables
		predictor_avg = paste0(predictor, "_avg")
		complementary_avg = paste0(predictors_dt[predictor, complementary], "_avg")

		## Colours
		colours = MetBrewer::met.brewer("Hiroshige", direction = ifelse(stri_detect(complementary_avg, regex = "^tas"), -1, 1))
		colours_str = gradColours(unlist(subset_trees[, ..complementary_avg]), colours, 30)

		# gradColours(seq(500, 1500, by = 100), colours, 30)

		# Plot 1: Scatter plot growth versus precip, with temperature as colour
		par(mar = c(5, 4, 2, 0))
		plot(x = unlist(subset_trees[, ..predictor_avg]), y = subset_trees[, annual_growth], ylim = c(0, 10),
			xlab = paste("Averaged", predictors_dt[predictor, explanation], "over the growing years"),
			ylab = "Growth in mm", pch = 20, col = paste0(colours_str, "66"))

		
		explanatory_var = unlist(subset_trees[, ..predictor_avg])
		m = mgcv::gam(subset_trees[, annual_growth] ~ s(explanatory_var))
		pred = cbind(x = sort(explanatory_var),
			as.data.frame(mgcv::predict.gam(m, se.fit = TRUE))[order(explanatory_var), ])
		lines(pred[, c("x", "fit")], col = "#9E391A", lwd = 2)
		polygon(c(pred$x, rev(pred$x)), c(pred$fit-1.96*pred$se.fit, rev(pred$fit+1.96*pred$se.fit)), col = "#354536", border = NA)

		# Plot 2: Top (precipitation) boxplot
		par(mar = c(0, 4, 0, 0))
		boxplot(subset_trees[, ..predictor_avg], xaxt = "n", yaxt = "n", bty = "n", col = "white", frame = FALSE, horizontal = TRUE)

		# Plot 3: Right (growth) boxplot
		par(mar = c(5, 0, 0, 0))
		boxplot(subset_trees[, annual_growth], xaxt = "n", yaxt = "n", bty = "n", col = "white", frame = FALSE)

		# Plot 4: Middle right (temperature) matrix gradient
		par(mar = c(5, 2, 2, 0))

		comp_min = round(min(floor(subset_trees[, ..complementary_avg]))) # Complementary min
		comp_max = round(max(ceiling(subset_trees[, ..complementary_avg]))) # Complementary max

		legend_image = grDevices::as.raster(matrix(gradColours(30:1, colours, 30), ncol = 1)) # 30:1 because raster starts from the top!
		plot(c(0, 2), c(comp_min, comp_max), type = "n", axes = FALSE, xlab = "", ylab = "",
			main = ifelse(predictors_dt[predictor, complementary] == "tas", "Avg temperature", "Avg precipitations"))
		if (predictors_dt[predictor, complementary] == "tas")
		{
			text(x = 1, y = seq(comp_min, comp_max, by = 3), labels = seq(comp_min, comp_max, by = 3))
		} else {
			text(x = 1, y = seq(comp_min, comp_max, by = 200), labels = seq(comp_min, comp_max, by = 200))
		}
		rasterImage(legend_image, 0, comp_min, 0.5, comp_max) # xleft, ybottom, xright, ytop

		dev.off()
	}
	print(paste(species, "done"))
}


#* -------------------------------------------------------------------------------------------
#! *******************************************************************************************
#? +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================            WHAT FOLLOWS ARE OLD STUFF            ================ ODOT#
#? +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! *******************************************************************************************
#* -------------------------------------------------------------------------------------------

climate_folder = ""
shapefile_folder = "/home/amael/shapefiles/"
epsg4326 = "EPSG:4326"

ls_folders = dir(treeData_folder)
countries_dt = data.table(country = c("France", "Deutschland"), trees_folder = c("FR IFN", ""),
	simplify = c(0.01, 100), shapefile = c("FRA_adm0", "germany"),
	epsg = c("2154", ""))

if (!all(countries_dt[, trees_folder] %in% ls_folders))
{
	warning(paste0("Folder not found for the following countries:\n- ",
		paste0(countries_dt[!(trees_folder %in% ls_folders), country], collapse = "\n- ")))
	
	countries_dt = countries_dt[trees_folder %in% ls_folders]
	warning("These countries were removed")
}

## Geographic zones
ls_countries = vector(mode = "list", length = countries_dt[, .N])
for (i in 1:countries_dt[, .N])
{
	country_shp = vect(paste0(shapefile_folder, countries_dt[i, country], "/", countries_dt[i, shapefile], ".shp"))
	ls_countries[[i]] = simplifyGeom(x = country_shp, tolerance = countries_dt[i, simplify])

	crs_country = setDT(crs(country_shp, describe = TRUE))
	if (paste0(crs_country[, authority], ":", crs_country[, code]) != epsg4326)
		ls_countries[[i]] = project(ls_countries[[i]], epsg4326)
}
names(ls_countries) = countries_dt[, country]

countries_dt[, nb_plots := 0]

countries_dt[, paste0(c("min", "quant_10", "median", "mean", "quant_90", "max"), "_T") := 0]
countries_dt[, paste0(c("min", "quant_10", "median", "mean", "quant_90", "max"), "_P") := 0]

temperature = rast(paste0(paste0(climate_folder, "tas/"), dir(path = paste0(climate_folder, "tas"), pattern = ".tif$")))
crs_temperature = setDT(crs(temperature, describe = TRUE))
if (paste0(crs_temperature[, authority], ":", crs_temperature[, code]) != epsg4326)
	temperature = project(temperature, epsg4326)

precipitation = rast(paste0(paste0(climate_folder, "pr/"), dir(path = paste0(climate_folder, "pr"), pattern = ".tif$")))
crs_precipitation = setDT(crs(precipitation, describe = TRUE))
if (paste0(crs_precipitation[, authority], ":", crs_precipitation[, code]) != epsg4326)
	precipitation = project(precipitation, epsg4326)

#### Plot data
i = 1
for (i in 1:countries_dt[, .N])
{
	## Load data
	treeData = readRDS(paste0(treeData_folder, "/trees_forest_reshaped.rds"))
	coords = unique(treeData[, .(xLambert93, yLambert93)])
	coords = vect(x = coords, geom = c("xLambert93", "yLambert93"), crs = paste0("EPSG:", countries_dt[i, epsg]))
	coords = project(x = coords, y = ls_countries[[i]])

	countries_dt[i, nb_plots := length(coords)]

	sub_temperature = setDT(extract(x = temperature, y = coords))
	sub_precipitation = setDT(extract(x = precipitation, y = coords))
	if ("ID" %in% names(sub_temperature))
	{
		warning("The column ID has been removed from the temperature data table")
		sub_temperature[, ID := NULL]
	}

	if ("ID" %in% names(sub_precipitation))
	{
		warning("The column ID has been removed from the precipitation data table")
		sub_precipitation[, ID := NULL]
	}

	countries_dt[i, paste0(c("min", "quant_10", "median", "quant_90", "max"), "_T") :=
		as.list(quantile(as.numeric(unlist(x = sub_temperature)), probs = c(0, 0.1, 0.5, 0.9, 1)))]
	countries_dt[i, paste0(c("min", "quant_10", "median", "quant_90", "max"), "_P") :=
		as.list(quantile(as.numeric(unlist(x = sub_precipitation)), probs = c(0, 0.1, 0.5, 0.9, 1)))]
	
	countries_dt[i, mean_T := mean(rowMeans(sub_temperature))]
	countries_dt[i, mean_P := mean(rowMeans(sub_precipitation))]
	
	pdf(paste0("./", countries_dt[i, country], ".pdf"))
	plot(ls_countries[[i]], axes = FALSE)
	points(coords, pch = 19, cex = 0.05)
	dev.off()
}

temperature_axis = seq(-100, 100, by = 5)
temperature_axis_min = temperature_axis[which.min(abs(temperature_axis - countries_dt[, min(min_T)]))]
temperature_axis_max = temperature_axis[which.min(abs(temperature_axis - countries_dt[, max(max_T)]))]
temperature_axis = seq(temperature_axis_min, temperature_axis_max, by = 5)

pdf("temperatureRange.pdf")
plot(0, pch = "", xlim = c(temperature_axis_min, temperature_axis_max), ylim = c(1, countries_dt[, .N]), axes = FALSE, bg = "transparent",
	xlab = "Temperature", ylab = "", asp = NA, type = "n")
for (i in 1:countries_dt[, .N])
{
	segments(x0 = countries_dt[i, min_T], y0 = i, x1 = countries_dt[i, max_T], y1 = i, col = "#21205B", lwd = 2)
	segments(x0 = countries_dt[i, quant_10_T], y0 = i, x1 = countries_dt[i, quant_90_T], y1 = i, col = "#D97A00", lwd = 6)
}

axis(side = 1, at = temperature_axis)
axis(side = 2, at = 1:countries_dt[, .N], labels = countries_dt[, country], las = 1, tick = FALSE)

dev.off()


precipitation_axis = seq(0, 10000, by = 5)
precipitation_axis_min = precipitation_axis[which.min(abs(precipitation_axis - countries_dt[, min(min_P)]))]
if (countries_dt[, min(min_P)] - precipitation_axis_min < 0)
	precipitation_axis_min = precipitation_axis_min - 5

precipitation_axis_max = precipitation_axis[which.min(abs(precipitation_axis - countries_dt[, max(max_P)]))]
if (countries_dt[, max(max_P)] - precipitation_axis_max > 0)
	precipitation_axis_max = precipitation_axis_max + 5

precipitation_axis = seq(precipitation_axis_min, precipitation_axis_max, by = 500)

pdf("precipitationRange.pdf")
plot(0, pch = "", xlim = c(precipitation_axis_min, precipitation_axis_max), ylim = c(1, countries_dt[, .N]), axes = FALSE,
	bg = "transparent", xlab = "Temperature", ylab = "", asp = NA, type = "n")
for (i in 1:countries_dt[, .N])
{
	segments(x0 = countries_dt[i, min_P], y0 = i, x1 = countries_dt[i, max_P], y1 = i, col = "#21205B", lwd = 2)
	segments(x0 = countries_dt[i, quant_10_P], y0 = i, x1 = countries_dt[i, quant_90_P], y1 = i, col = "#D97A00", lwd = 6)
}

axis(side = 1, at = precipitation_axis)
axis(side = 2, at = 1:countries_dt[, .N], labels = countries_dt[, country], las = 1, tick = FALSE)

dev.off()
