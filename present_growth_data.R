
#### Aim of prog: Describe the growth data
## List of plots:
#	1. Climate range for few species (most abundant)
#	2. Number of individuals in the database per species and country
#	3. Maps of plots
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

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
			if (is.list(temporary) & all(names(temporary) == c("time_space", "env", "loaded")))
			{
				env = temporary[["env"]]
			} else {
				collectEnvironment(trees_dt, env_dt, climVars, file = NULL) # Recursive call with option NULL to force computation
				warning("The loaded object did not have the required format. Recomputing")
			}
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

computeEnvironment = function(trees_dt, env, ls_species = NULL, climVars = c("pr", "tas"))
{
	if (!is.null(ls_species))
		trees_dt = trees_dt[speciesName_sci %in% ls_species]

	ls_species = trees_dt[, unique(speciesName_sci)]

	quantiles_ls = vector(mode = "list", length = length(ls_species))
	names(quantiles_ls) = ls_species
	setkey(env, plot_id, year)

	for (species in ls_species)
	{
		time_space = unique(trees_dt[speciesName_sci == species, .(plot_id, year)])
		time_space[, year_start := min(year), by = plot_id]
		time_space[, year_end := max(year), by = plot_id]
		time_space[, year := NULL]
		time_space = unique(time_space)

		sp_clim = vector(mode = "list", length = time_space[, .N])

		for (i in 1:time_space[, .N])
		{
			current_plot = time_space[i, plot_id]
			year_start = time_space[i, year_start]
			year_end = time_space[i, year_end]
			sp_clim[[i]] = env[.(current_plot, year_start:year_end)]
		}
		sp_clim = rbindlist(sp_clim)

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

climBounds = computeEnvironment(treeData, env[["env"]],
	ls_species = c("Betula pendula", "Fagus sylvatica", "Pinus sylvestris", "Quercus robur"))

