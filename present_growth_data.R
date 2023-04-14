
#### Aim of prog: Describe the growth data
## List of plots:
#	1. Climate range and number of individuals in the data set per species and country
#	2. Maps of plots
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)
library(terra)

if (!nzchar(system.file(package = "MetBrewer")))
	install.packages("MetBrewer")

#? --------------------------------------------------------------------------------------------------------
######## PART 1: Climate range and number of individuals in the data set per species and country
#? --------------------------------------------------------------------------------------------------------

#### Tool functions
collectEnvironment = function(trees_dt, env_dt, climVars = c("pr", "tas"), file = NULL)
{
	time_space = unique(trees_dt[, .(plot_id, year)])
	time_space[, year_start := min(year), by = plot_id]
	time_space[, year_end := max(year), by = plot_id]
	time_space[, year := NULL]
	time_space = unique(time_space)

	setkey(env_dt, plot_id, year)
	
	loaded = FALSE
	env = vector(mode = "list", length = time_space[, .N])
	names(env) = time_space[, paste(plot_id, year_start, year_end, sep = "_")]
	
	if (!is.null(file))
	{
		if (file.exists(file))
		{
			temporary = readRDS(file)
			loaded = TRUE
			warning(paste("Found file", file, "and loaded it rather than computing the environment"))
			if (is.list(temporary) && all(names(temporary) == c("time_space", "env", "loaded")))
			{
				env = temporary[["env"]]
			} else {
				collectEnvironment(trees_dt, env_dt, climVars, file = NULL) # Recursive call with option NULL to force computation
				warning("The loaded object did not have the required format. Recomputing")
			}
		} else {
			envList = collectEnvironment(trees_dt, env_dt, climVars, file = NULL) # Recursive call with option NULL to force computation
			return (list(time_space = time_space, env = envList[["env"]], loaded = loaded))
		}
	} else {
		for (i in 1:time_space[, .N])
		{
			current_plot = time_space[i, plot_id]
			year_start = time_space[i, year_start]
			year_end = time_space[i, year_end]
			env[[i]] = env_dt[.(current_plot, year_start:year_end)]
			if (i %% 500 == 0)
				print(paste0(round(i*100/time_space[, .N], 1), "% done"))
		}
		print("100% done")
		env = rbindlist(env, idcol = "plot_id_startYear_endYear")
	}
	return (list(time_space = time_space, env = env, loaded = loaded))
}

computeEnvironment = function(trees_dt, env, ls_species = NULL, climVars = c("pr", "tas"), spClim_folder,
	loading = TRUE, saving = TRUE, overwriting = FALSE)
{
	if (!is.null(ls_species))
		trees_dt = trees_dt[speciesName_sci %in% ls_species]

	ls_species = trees_dt[, unique(speciesName_sci)]

	quantiles_ls = vector(mode = "list", length = length(ls_species))
	names(quantiles_ls) = ls_species
	setkey(env, plot_id, year)

	if (!dir.exists(spClim_folder))
	{
		if (loading)
		{
			warning("Loading is set to TRUE, but the folder does not exist! loading set to FALSE, folder created and saving set to TRUE")
			loading = FALSE
			saving = TRUE
		}
		dir.create(spClim_folder)
	}

	for (species in ls_species)
	{
		time_space = unique(trees_dt[speciesName_sci == species, .(plot_id, year)])
		time_space[, year_start := min(year), by = plot_id]
		time_space[, year_end := max(year), by = plot_id]
		time_space[, year := NULL]
		time_space = unique(time_space)

		sp_clim_file = paste0(spClim_folder, species, ".rds")

		if (loading && file.exists(sp_clim_file))
		{
			sp_clim = readRDS(sp_clim_file)
		} else {
			sp_clim = vector(mode = "list", length = time_space[, .N])

			for (i in 1:time_space[, .N])
			{
				current_plot = time_space[i, plot_id]
				year_start = time_space[i, year_start]
				year_end = time_space[i, year_end]
				sp_clim[[i]] = env[.(current_plot, year_start:year_end)]
			}
			sp_clim = rbindlist(sp_clim)
		}

		if (saving)
		{
			if (file.exists(sp_clim_file))
			{				
				if (overwriting)
				{
					saveRDS(sp_clim, sp_clim_file)
				} else {
					if (!loading)
						warning(paste0("The file <", sp_clim_file, "> already exists. To overwrite it, use the option overwriting = TRUE"))
				}
			} else {
				saveRDS(sp_clim, sp_clim_file)
				print(paste0("File saved at <", sp_clim_file, ">"))
			}
		}

		quantiles_ls[[species]] = sp_clim[, lapply(.SD, function (x) {c(quantile(x, c(0, 0.05, 0.5, 0.95, 1)), mean(x))}),
			.SDcols = climVars]
		quantiles_ls[[species]][, stats := c("min", "p05", "med", "p95", "max", "avg")]
		print(paste("Species", species, "done"))
	}
	quantiles_ls = rbindlist(quantiles_ls, idcol = "speciesName_sci")
	return(quantiles_ls)
}

#### Load data
## Common variables
dataFolder = "/home/amael/project_ssm/inventories/growth/"
spClim_folder = "/home/amael/project_ssm/inventories/growth/spClim/"

## Tree data
treeData = readRDS(paste0(dataFolder, "standardised_european_growth_data_reshaped.rds"))

## Climate data
climate = readRDS(paste0(dataFolder, "europe_reshaped_climate.rds"))

## Soil data (pH) data
soil = readRDS(paste0(dataFolder, "europe_reshaped_soil.rds"))

## Interpolated basal area data
standBasalArea = readRDS(paste0(dataFolder, "europe_reshaped_standBasalArea.rds"))

#### Info species
## Number of individuals per species, country
species_data = treeData[, length(unique(tree_id)), by = .(speciesName_sci, plot_id, country)]
species_data = species_data[, sum(V1), by = .(speciesName_sci, country)]
setnames(species_data, old = "V1", new = "sp_indivs_per_country")
species_data[, sp_indivs_tot := sum(sp_indivs_per_country), by = speciesName_sci]

## Climate data per species
env = collectEnvironment(treeData, climate, climVars = c("pr", "tas"), file = paste0(dataFolder, "time_space.rds"))

if (!env[["loaded"]])
	saveRDS(env, paste0(dataFolder, "time_space.rds"))

climBounds = computeEnvironment(treeData, env[["env"]], spClim_folder = spClim_folder)

setkey(climBounds, speciesName_sci, stats)
saveRDS(climBounds, paste0(dataFolder, "climBounds.rds"))

#### Plot climate data per species
## Common variables
keptSpecies = c("Abies alba", "Acer campestre", "Acer platanoides", "Acer pseudoplatanus", "Alnus glutinosa", "Alnus incana",
		"Betula pendula", "Betula pubescens", "Fagus sylvatica", "Larix decidua", "Larix kaempferi", "Picea abies", "Picea sitchensis",
		"Pinus contorta", "Pinus halepensis", "Pinus nigra", "Pinus pinaster", "Pinus sylvestris", "Populus tremula", "Pseudotsuga menziesii",
		"Quercus petraea", "Quercus robur", "Salix caprea", "Sorbus aucuparia", "Tilia cordata", "Ulmus minor")

climBounds = climBounds[.(keptSpecies)]
ls_species = climBounds[, unique(speciesName_sci)] # i.e., keptSpecies

ls_vars = colnames(climBounds)[!(colnames(climBounds) %in% c("speciesName_sci", "stats"))]
y_max = length(ls_species) - 1

for (currentVar in ls_vars)
{
	# pdf(paste0(currentVar, ".pdf"), height = 6, width = 9.708204)
	tikz(file = paste0(currentVar, ".tex"), width = 5.562306, height = 9) # Width = 18 cm or 7.08 Inches, Golden ratio
	x_min = min(climBounds[stats == "p05", ..currentVar])
	x_max = max(climBounds[stats == "p95", ..currentVar])

	counter = y_max

	plot(0, type = "n", axes = FALSE, ann = FALSE, xlim = c(x_min, x_max), ylim = c(0, y_max + 0.5))
	par(mar = c(3, 3, 2, 2), oma = c(0, 8, 0, 3))

	for (species in ls_species)
	{
		segments(x0 = unlist(climBounds[.(species, "p05"), ..currentVar]), y0 = counter,
			x1 = unlist(climBounds[.(species, "p95"), ..currentVar]))
		segments(x0 = unlist(climBounds[.(species, "p05"), ..currentVar]), y0 = counter - 0.05, y1 = counter + 0.05, lwd = 2)
		points(x = unlist(climBounds[.(species, "med"), ..currentVar]), y = counter, col = "#95C36E", bg = "#95C36E", pch = 23)
		segments(x0 = unlist(climBounds[.(species, "p95"), ..currentVar]), y0 = counter - 0.05, y1 = counter + 0.05, lwd = 2)
		points(x = unlist(climBounds[.(species, "avg"), ..currentVar]), y = counter, col = "#295384", pch = 19)

		counter = counter - 1
	}
	axis(side = 1, at = round(seq(x_min, x_max, length.out = 3), 1),
		labels = round(seq(x_min, x_max, length.out = 3), ifelse(currentVar == "pr", 0, 1)))
	axis(side = 2, at = y_max:0, labels = ls_species, las = 1)
	legend(x = "topleft", xpd = TRUE, legend = c("Median", "Mean"), pch = c(23, 19), col = c("#95C36E", "#295384"),
		pt.bg = c("#95C36E", "#295384"), bty = "n", horiz = TRUE, inset = c(0, -0.01))
	dev.off()
}

#### Table climate range
climBounds_xt = climBounds[.(unique(speciesName_sci), rep(c("min", "max"), each = length(ls_species)))]
climBounds_xt = dcast(climBounds_xt, speciesName_sci ~ stats, value.var = c("pr", "tas"))
climBounds_xt = merge(x = climBounds_xt, y = unique(species_data[, .(speciesName_sci, sp_indivs_tot)]), by = "speciesName_sci")

setcolorder(x = climBounds_xt, neworder = c("speciesName_sci", "sp_indivs_tot", paste(rep(ls_vars, each = 2), c("min", "max"), sep = "_")))

climBounds_xt = xtable(x = climBounds_xt, align = c("r", "r", "c", "c", "c", "c", "c"), digits = 1)
print(x = climBounds_xt, file = "./climBounds.tex", include.rownames = FALSE, table.placement = "h", booktabs = TRUE)

#? --------------------------------------------------------------------------------------------------------
######## PART 3: Map of plots
#? --------------------------------------------------------------------------------------------------------
#### Load data
## Folders
shapefile_folder = "/home/amael/shapefiles/europe/continent/"

## Tree coords and shapefile
# Get data
coords = vect(unique(treeData[, .(x, y)]), geom = c("x", "y"), crs = "EPSG:4326")
europe = vect(paste0(shapefile_folder, "europe.shp"))

# Bounding box
bbox_europe = c(xmin = -5.69, xmax = 28, ymin = 38, ymax = 69.53)
europe_extent = ext(bbox_europe)

# Projecting to EPSG:3035 (Lambert azimuthal equal-area projection), recommended by European Environment Agency
coords = project(coords, "EPSG:3035")
europe = project(europe, "EPSG:3035")
europe_extent = project(rast(europe_extent, crs = "EPSG:4326"), "EPSG:3035")

europe = crop(x = europe, y = europe_extent)

#### Plot
# pdf("plots_location_growth_subsample_2.pdf", height = 10, width = 10)
pdf("plots_location_growth.pdf", height = 10, width = 10)
plot(europe, col = "#C4AC7C44", border = "#9E391A", axes = FALSE)
plot(coords, pch = 20, cex = 0.025, col = "#354536", add = TRUE)
plot(europe, col = NA, border = "#9E391A", add = TRUE)
dev.off()

#### Number of plots per country
nbPlotsPerCountry = unique(treeData[, .(x, y, country)])[, .N, by = country]
setkey(nbPlotsPerCountry, country)

print(paste("Total number of plots:", nbPlotsPerCountry[, sum(N)]))
print(paste("Total number of individuals:", species_data[, sum(sp_indivs_per_country)]))
print(paste("Total number of species:", length(keptSpecies)))

nbPlotsPerCountry_xt = xtable(x = nbPlotsPerCountry, align = c("r", "r", "c"), digits = 1)
print(x = nbPlotsPerCountry_xt, file = "./nbPlotsPerCountry.tex", include.rownames = FALSE, table.placement = "h", booktabs = TRUE)
