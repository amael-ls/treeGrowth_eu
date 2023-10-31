
#### Aim of prog: Show the whole distribution of trees, how much is covered by our dataset (spatial and climatic extent)
## Comments:
# The chorological maps have been downloaded from:
#	Caudullo, Giovanni and Erik, Welk and Jesús, San-Miguel-Ayanz (2017)
#	“Chorological Maps for the Main European Woody Species.” in Data in Brief 12: 662–66.
#	https://doi.org/10.1016/j.dib.2017.05.007.
#
# The european shapefile is from:
#	https://ec.europa.eu/eurostat/web/gisco/
#
# The world shapefile is from:
#	https://datacatalog.worldbank.org/search/dataset/0038272/World-Bank-Official-Boundaries
#
# The ocean shapefile is from:
#	https://www.naturalearthdata.com/downloads/50m-physical-vectors/50m-ocean/

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(stringi)
library(terra)

#### Compute climate range for each species
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
ls_var = c("tas", "pr")

dir_shapefiles = "/home/amael/shapefiles/trees/europe/chorological_maps/"

europe = vect("/home/amael/shapefiles/europe/continent/europe.shp")
world = vect("/home/amael/shapefiles/world/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp")
ocean = vect("/home/amael/shapefiles/world/ne_50m_ocean/ne_50m_ocean.shp")

europe = project(europe, world)
ocean = project(ocean, world)

ls_countries = c(europe[europe$NAME_LATN == "France"],
	europe[europe$NAME_LATN == "Deutschland"],
	europe[europe$NAME_LATN == "Sverige"])
freq_hatch = data.table(speciesName_sci = ls_species, n = c(3, 8, 8, 12, 3, 8), height = c(4.5, 4.5, 4.5, 4.5, 4.5, 4.5),
	width = c(9, 4.5, 4.5, 4.5, 9, 6), key = "speciesName_sci") # Empirical values for n

clim_range = vector(mode = "list", length = length(ls_species))
names(clim_range) = ls_species

areaCover = readRDS("./areaCover.rds")

## Load species shapefiles and world climate rasters to compute range
for (species in ls_species)
{
	species_underscore = stri_replace(str = species, replacement = "_", regex = " ")
	currentDir = paste0(dir_shapefiles, species, "/shapefiles/")
	shapefile_name = list.files(path = currentDir, pattern = paste0("^", species_underscore, "_plg.shp"))

	shp = project(vect(paste0(currentDir, shapefile_name)), world)
	if (!all(is.valid(shp)))
		shp = makeValid(shp)
	
	shp = crop(x = shp, y = world)

	clim_range[[species]] = data.table(variable = ls_var, min = -Inf, q025 = Inf, q50 = Inf, q975 = Inf, max = Inf, key = "variable")
	for (currentVariable in ls_var)
		clim_range[[species]][currentVariable, c("min", "q025", "q50", "q975", "max") :=
			as.list(readRDS(paste0("./", species, "/", currentVariable, "_qt.rds"))[c("0%", "2.5%", "50%", "97.5%", "100%")])]

	hatch = vector(mode = "list", length = length(ls_countries))

	for (country_id in seq_along(ls_countries))
	{
		bbox = ext(ls_countries[country_id])

		span_x = bbox$xmax - bbox$xmin
		span_y = bbox$ymax - bbox$ymin

		delta_x = span_x/freq_hatch[species, n]
		delta_y = span_y/freq_hatch[species, n]

		ls_coords = matrix(data = NA, nrow = 4*freq_hatch[species, n], ncol = 4)
		colnames(ls_coords) = c("object", "part", "x", "y")

		for (i in seq_len(freq_hatch[species, n]))
		{
			ls_coords[4*i - 3, ] = c(2*i, 1, bbox$xmin + i*delta_x, bbox$ymin)
			ls_coords[4*i - 2, ] = c(2*i, 1, bbox$xmin, bbox$ymin + i*delta_y)
			ls_coords[4*i - 1, ] = c(2*i + 1, 1, bbox$xmin + i*delta_x, bbox$ymax)
			ls_coords[4*i, ] = c(2*i + 1, 1, bbox$xmax, bbox$ymin + i*delta_y)
		}

		hatch[[country_id]] = vect(ls_coords, "lines", crs = crs(world))
		hatch[[country_id]] = crop(x = hatch[[country_id]], y = ls_countries[country_id])
		hatch[[country_id]] = crop(x = hatch[[country_id]], y = shp)
	}

	print(paste0("Overlap ", species, ": ", areaCover[species, round(percentCover, 2)], "%"))
	pdf(paste0("./distribution_", species, ".pdf"), height = freq_hatch[species, height], width = freq_hatch[species, width])
	plot(shp, col = "#3F7752AA", main = "", axes = FALSE, mar = c(0, 0, 0, 0))
	plot(world, add = TRUE)
	plot(ocean, add = TRUE, col = "#0C93A344")
	for (i in seq_along(hatch))
		if (!all(dim(hatch[[i]]) == 0))
			plot(hatch[[i]], add = TRUE, lwd = 1.8, col = "red")
	dev.off()
}

#### Plot climate range of the distribution (clim_range) and of the data (climBounds)
clim_range = rbindlist(clim_range, idcol = "speciesName_sci")
setkey(clim_range, speciesName_sci, variable)

climBounds = readRDS("/home/amael/project_ssm/inventories/growth/climBounds.rds")
climBounds = climBounds[.(ls_species)]
y_max = length(ls_species) - 1

delta = 0.1

for (currentVariable in ls_var)
{
	# pdf(paste0(currentVariable, ".pdf"), height = 6, width = 9.708204)
	tikz(file = paste0(currentVariable, ".tex"), width = 7, height = 9) # Width = 18 cm or 7.08 Inches, Golden ratio
	x_min = min(climBounds[stats == "p05", ..currentVariable], clim_range[.(unique(speciesName_sci), currentVariable), q025])
	x_max = max(climBounds[stats == "p95", ..currentVariable], clim_range[.(unique(speciesName_sci), currentVariable), q975])

	counter = y_max

	plot(0, type = "n", axes = FALSE, ann = FALSE, xlim = c(x_min, x_max), ylim = c(0, y_max + 0.5))
	par(mar = c(3, 3, 2, 4), oma = c(0, 8, 0, 3))

	for (species in ls_species)
	{
		# From plot data
		segments(x0 = unlist(climBounds[.(species, "p05"), ..currentVariable]), y0 = counter,
			x1 = unlist(climBounds[.(species, "p95"), ..currentVariable]))
		segments(x0 = unlist(climBounds[.(species, "p05"), ..currentVariable]), y0 = counter - 0.05, y1 = counter + 0.05, lwd = 2)
		points(x = unlist(climBounds[.(species, "med"), ..currentVariable]), y = counter, col = "#295384", bg = "#295384", pch = 23)
		segments(x0 = unlist(climBounds[.(species, "p95"), ..currentVariable]), y0 = counter - 0.05, y1 = counter + 0.05, lwd = 2)
		
		# From distribution
		segments(x0 = unlist(clim_range[.(species, currentVariable), q025]), y0 = counter + delta,
			x1 = unlist(clim_range[.(species, currentVariable), q975]), lty = 2)
		segments(x0 = unlist(clim_range[.(species, currentVariable), q025]), y0 = counter - 0.05 + delta, y1 = counter + 0.05 + delta, lwd = 2)
		points(x = unlist(clim_range[.(species, currentVariable), q50]), y = counter + delta, col = "#295384", bg = "#295384", pch = 23)
		segments(x0 = unlist(clim_range[.(species, currentVariable), q975]), y0 = counter - 0.05 + delta, y1 = counter + 0.05 + delta, lwd = 2)

		counter = counter - 1
	}
	axis(side = 1, at = round(seq(x_min, x_max, length.out = 3), 1),
		labels = round(seq(x_min, x_max, length.out = 3), ifelse(currentVariable == "pr", 0, 1)))
	axis(side = 2, at = y_max:0, labels = ls_species, las = 1)
	legend(x = "topleft", legend = c("Median", "Distribution", "Data"), lty = c(NA, 2, 1), pch = c(23, NA, NA), xpd = TRUE,
		col = c("#295384", "black", "black"), pt.bg = c("#295384", NA, NA), bty = "n", horiz = TRUE, inset = c(0, -0.01))
	dev.off()
}

#? --------------------------------------------------------------------------------------------------------
######## PART II: Spatial plot of the data (plot location)
#? --------------------------------------------------------------------------------------------------------
#### Load data
## Folders
treeData_folder = "/home/amael/project_ssm/inventories/growth/"

## Tree data
treeData = readRDS(paste0(treeData_folder, "standardised_european_growth_data_reshaped.rds"))
coords = vect(unique(treeData[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:4326")

# coords = sample(coords, 45000) # Useful for presentations where a too dense picture shows nothing!

#### Plot
# pdf("plots_location_growth_subsample_2.pdf", height = 10, width = 10)
pdf("plots_location_growth.pdf", height = 10, width = 10)
plot(europe, col = "#C4AC7C44", border = "#9E391A", axes = FALSE)
plot(coords, pch = 20, cex = 0.025, col = "#354536", add = TRUE)
plot(europe, col = NA, border = "#9E391A", add = TRUE)
dev.off()