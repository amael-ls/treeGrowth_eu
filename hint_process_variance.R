
#### Aim of prog: Get a hint on the variance of growth
## Comments:
# I check in this program what is the variance of growth for individuals that are similar, and in similar conditions.
# I measure the similarity with the Euclidean distance between two individuals in the space Environment x dbh x species

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)

#### Tool functions
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
	dt[, (col) := (shift(dbh, n = 1, type = "lead", fill = NA) - dbh)/(shift(year, n = 1, type = "lead", fill = NA) - year), by = byCols]
	
	if (!("deltaYear" %in% names(dt)))
		dt[, deltaYear := shift(year, n = 1, type = "lead", fill = NA) - year, by = byCols]
}

## Function to aggregate environment variables over a period and relate them to individuals
computeEnvironment = function(trees_dt, env_dt, climVars = c("pr", "tas"))
{
	time_space = unique(trees_dt[, .(plot_id, year, deltaYear)])
	time_space[, year_end := year + deltaYear]
	time_space[, deltaYear := NULL]
	
	varsComputation = c("min", "max", "mean", "sd")
	newCols = CJ(climVars, varsComputation, sorted = FALSE)[, paste(climVars, varsComputation, sep ="_")]
	
	values_ls = vector(mode = "list", length = 4)
	names(values_ls) = c("min", "max", "mean", "sd")
	
	time_space[, c(newCols) := 0]
	
	setnames(time_space, old = "year", new = "year_start")
	
	if (file.exists("timeSpace_germany.rds"))
	{
		time_space = readRDS("timeSpace_germany.rds")
		warning("Found file `timeSpace_germany.rds` and loaded it rather than computing the environment")
	} else {
		for (i in 1:time_space[, .N])
		{
			current_plot = time_space[i, plot_id]
			year_start = time_space[i, year_start]
			year_end = time_space[i, year_end]
			temporary = env_dt[(plot_id == current_plot) & (year_start <= year) & (year <= year_end)]
			
			values_ls[["min"]] = temporary[, lapply(.SD, min), .SDcols = climVars]
			values_ls[["max"]] = temporary[, lapply(.SD, max), .SDcols = climVars]
			values_ls[["mean"]] = temporary[, lapply(.SD, mean), .SDcols = climVars]
			values_ls[["sd"]] = temporary[, lapply(.SD, sd), .SDcols = climVars]
			values = rbindlist(values_ls)
			
			time_space[i, c(newCols) := as.list(unlist(values))]
			
			if (i %% 500 == 0)
				print(paste0(round(i*100/time_space[, .N], 1), "% done"))
		}
		print("100% done")
	}
	return (list(time_space = time_space, varnames = newCols))
}

#### Load data
## Paths
mainFolder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

clim_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

soil_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(soil_folder))
	stop(paste0("Folder\n\t", soil_folder, "\ndoes not exist"))

standBasalArea_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(standBasalArea_folder))
	stop(paste0("Folder\n\t", standBasalArea_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))
treeData = treeData[country == "germany"] # Only keep Germany since I currently do not work with France (different error structure)

## Compute growth
treeData = computeDiametralGrowth(treeData, byCols = c("speciesName_sci", "plot_id", "tree_id")) # Ok because plot_id is country specific
growth_dt = na.omit(treeData)

## Read climate
climate = readRDS(paste0(clim_folder, "europe_reshaped_climate.rds"))

## Read soil data (pH)
soil = readRDS(paste0(soil_folder, "europe_reshaped_soil.rds"))

## Read interpolated basal area data
standBasalArea = readRDS(paste0(standBasalArea_folder, "europe_reshaped_standBasalArea.rds"))

#### Merge climate and growth data tables
## Compute environmental variables (range, mean, and sd per plot per individual timespan)
env = computeEnvironment(growth_dt, climate, climVars = c("pr", "tas"))
time_space = env[["time_space"]]
envCols = env[["varnames"]]
envCols_dbh = c(envCols, "dbh")

## Change names growth for merge
setnames(growth_dt, old = "year", new = "year_start")
growth_dt[, year_end := year_start + deltaYear]

## Merging
growth_dt = growth_dt[time_space, on = c("plot_id", "year_start", "year_end")]

#### Slicing
## Dbh slicing (based on a visual histogram)
slice_1 = growth_dt[(267 < dbh) & (dbh < 272)]

## Precipitation slicing (based on a visual histogram)
slice_1 = slice_1[(612 < pr_mean) & (pr_mean < 622)]

## Temperature (tas) slicing
hist(slice_1[, tas_mean]) # no need to slice, I condider it fine

slice_1[, mean(log(growth))]
slice_1[, sd(log(growth))]

hist(slice_1[, growth], probability = TRUE, ylim = c(0, 0.27))
curve(dlnorm(x, meanlog = slice_1[, mean(log(growth))], sdlog = slice_1[, sd(log(growth))]), add = TRUE, col = "orange", lwd = 2)

slice_1[, mean(growth)]
slice_1[, sd(growth)]

curve(dgamma(x, shape = slice_1[, mean(growth)]^2/slice_1[, var(growth)],
	rate = slice_1[, mean(growth)]/slice_1[, var(growth)]), add = TRUE, col = "red", lwd = 2)

slice_1[, .N, by = speciesName_sci]

## Pinus sylvestris
slice_1 = slice_1[speciesName_sci == "Pinus sylvestris"]

slice_1[, mean(log(growth))]
slice_1[, sd(log(growth))]

hist(slice_1[, growth], probability = TRUE, breaks = seq(0, 7, length.out = 15))
curve(dlnorm(x, meanlog = slice_1[, mean(log(growth))], sdlog = slice_1[, sd(log(growth))]), add = TRUE, col = "orange", lwd = 2)

slice_1[, mean(growth)]
slice_1[, sd(growth)]

curve(dgamma(x, shape = slice_1[, mean(growth)]^2/slice_1[, var(growth)],
	rate = slice_1[, mean(growth)]/slice_1[, var(growth)]), add = TRUE, col = "red", lwd = 2)

slice_1[, .N, by = speciesName_sci]

slice_1[, growth_scaled := growth/sd(dbh)]

hist(slice_1[, growth_scaled], probability = TRUE, breaks = seq(0, 7, length.out = 15))
curve(dlnorm(x, meanlog = slice_1[, mean(log(growth_scaled))], sdlog = slice_1[, sd(log(growth_scaled))]),
	add = TRUE, col = "orange", lwd = 2)

curve(dgamma(x, shape = slice_1[, mean(growth_scaled)]^2/slice_1[, var(growth_scaled)],
	rate = slice_1[, mean(growth_scaled)]/slice_1[, var(growth_scaled)]), add = TRUE, col = "red", lwd = 2)

#### Compute distances between individuals
## Transform all variables between 0 and 1. This transformation preserves the order of Euclidean distance (2-norm in R^k)
growth_dt[, c(envCols_dbh) := lapply(.SD, function(x) {(x - min(x))/(max(x) - min(x))}), .SDcols = envCols_dbh]

## Compute distance between individuals
# Remove trees measured more than twice (to avoid comparing two growth from the same individual)
growth_dt[, n_growth := .N, by = .(tree_id, plot_id)]
growth_dt = growth_dt[n_growth == 1]
growth_dt[, n_growth := NULL]

# Common variable
ls_species = sort(growth_dt[, unique(speciesName_sci)])
ls_species = c("Betula pubescens", "Fagus sylvatica", "Pinus sylvestris", "Quercus robur")

distance_ls = vector(mode = "list", length = length(ls_species))
names(distance_ls) = ls_species

# Distances (stored in a list of matrices distances)
for (species in ls_species)
	distance_ls[[species]] = dist(growth_dt[speciesName_sci == species, ..envCols_dbh], method = "euclidean", diag = FALSE,	
		upper = FALSE, p = 2)

species = ls_species[1]
speciesDistance = distance_ls[[species]]

#### Access data and plot distribution growth for similar individuals (similar dbh and environment)
## Compute number of individuals per species
individuals_per_species = growth_dt[, .N, by = speciesName_sci]
setkey(individuals_per_species, speciesName_sci)

if (individuals_per_species[, sum(N)] != growth_dt[, .N])
	warning("Mismatch between total number of individuals and species-specific numbers.")

n_indiv = individuals_per_species[species, N]
indiv1 = 1
indiv2 = 2:n_indiv

dissimilarity_1_to_others = speciesDistance[n_indiv*(indiv1 - 1) - indiv1*(indiv1 - 1)/2 + indiv2 - indiv1]
range(dissimilarity_1_to_others)

threshold = quantile(dissimilarity_1_to_others, seq(0, 1, 0.05))["5%"]

others = which(dissimilarity_1_to_others <= threshold) - n_indiv*(indiv1 - 1) + indiv1*(indiv1 - 1)/2 + indiv1

species_dt = growth_dt[speciesName_sci == species]
species_dt[c(indiv1, others)]


#### Crash test zone
sd_dbh = 150
curve(dgamma(x, shape = 1.5/1.0, rate = sd_dbh*sqrt(1.5)/1.0), to = 0.1)

curve(dgamma(x, shape = 1.5/1.0, rate = sqrt(1.5)/1.0), to = 5)
