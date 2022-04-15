#### Aim of prog: Plot tree data on a map and few others informations (climate range, etc)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#### Load data
## Common variables
treeData_folder = "/home/amael/project_ssm/inventories/growth/"
shapefile_folder = "/home/amael/shapefiles/europe/"

treeData = readRDS(paste0(treeData_folder, "standardised_european_growth_data_reshaped.rds"))
coords = vect(unique(treeData[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:4326")

europe = vect(paste0(shapefile_folder, "europe.shp"))

pdf("plots_location_growth.pdf", height = 10, width = 10)
# plot(0, pch = "", xlim = ext(europe)[1:2], ylim = ext(europe)[3:4], axes = FALSE, bg = "transparent",
	# xlab = "", ylab = "", asp = NA, type = "n")
plot(europe, col = "#C4AC7C44", border = "#9E391A", axes = FALSE)
plot(coords, pch = 19, cex = 0.025, col = "#354536", add = TRUE)
plot(europe, col = NA, border = "#9E391A", add = TRUE)
dev.off()






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
