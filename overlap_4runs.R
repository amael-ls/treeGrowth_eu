
#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#### Tool functions
## Function to compute all the combinations of overlaps for one parameter
overlap_fct = function(polygons_ls, n_runs)
{
	checkPolygons = sapply(polygons_ls, terra::is.valid)

	if (!all(checkPolygons))
	{
		nn = names(polygons_ls)
		warning(paste("The polygons", nn[!checkPolygons], "were not valid, but have been forced to be via terra::makeValid"))
		for (polyg in nn[!checkPolygons])
			polygons_ls[[polyg]] = terra::makeValid(polygons_ls[[polyg]])
	}

	overlap_ls = vector(mode = "list", n_runs - 1)
	names(overlap_ls) = paste0(2:n_runs, "_elements")

	for (nbElements in 2:n_runs) # Number of elements in the combination
	{
		combinations_ls = combn(x = 1:n_runs, m = nbElements, simplify = TRUE)
		results_dt = data.table(combination = character(ncol(combinations_ls)), overlap = numeric(ncol(combinations_ls)))

		count = 1
		
		for (j in seq_len(ncol(combinations_ls))) # Loop among the combinations with 'nbElements' runs
		{
			selected_runs = paste0("run_", combinations_ls[, j])
			currentCombination = paste(combinations_ls[, j], collapse = "_")

			union_polygon = polygons_ls[[selected_runs[1]]]
			inter_polygon = polygons_ls[[selected_runs[1]]]

			for (currentRun in selected_runs[2:length(selected_runs)]) # Loop among the selected runs
			{
				union_polygon = terra::union(union_polygon, polygons_ls[[currentRun]])
				inter_polygon = terra::intersect(inter_polygon, polygons_ls[[currentRun]])
			}
			union_polygon = terra::aggregate(union_polygon)
			results_dt[count, combination := currentCombination]
			suppressWarnings(results_dt[count, overlap := terra::expanse(inter_polygon)/terra::expanse(union_polygon)])
			count = count + 1
		}
		overlap_ls[[paste0(nbElements, "_elements")]] = results_dt
	}


	overlap_dt = rbindlist(overlap_ls, idcol = "elements")
	overlap_dt[, elements := stri_sub(str = elements, to = stri_locate_first(str = elements, regex = "_")[, "start"] - 1)]
	
	return(overlap_dt)
}

## Other functions
source("./toolFunctions.R")

#### Common variables
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
	stop("Supply the species", call. = FALSE)

species = as.character(args[1])
path = paste0("./", species, "/")

n_runs = 4 # Number of runs used in growth_subsample.R

## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

threshold_indiv = 12000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2023/01/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_de-fr-sw_12000_main.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

nb_nfi = infoSpecies[species, n_nfi]

## Parameters
params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "sigmaProc", "etaObs")

if (nb_nfi > 1)
	params = expand(params, nb_nfi)[["new_names"]]

n_params = length(params)

## Posteriors object
polygons_ls = vector(mode = "list", length = n_params)
names(polygons_ls) = params

for (currentParam in params)
{
	polygons_ls[[currentParam]] = vector(mode = "list", length = n_runs)
	names(polygons_ls[[currentParam]]) = paste0("run_", 1:n_runs)
}

#### Compute overlap posteriors
## Load results and extract posteriors
for (i in 1:n_runs)
{
	info_lastRun = getLastRun(path = path, extension = "_main.rds$", run = i)
	lastRun = info_lastRun[["file"]]
	results = readRDS(paste0(path, lastRun))

	currentRun = paste0("run_", i)

	for (currentParam in params)
	{
		currentDensity = density(results$draws(currentParam), n = 2048)
		x = currentDensity$x
		y = currentDensity$y

		temporary = cbind(id = 1, part = 1, x, y)
		polygons_ls[[currentParam]][[currentRun]] = vect(temporary, type = "polygons")
	}
	
	print(paste("Run", i, "done"))
}

## Overlaps for all the combinations
overlap_ls = vector(mode = "list", length(params))
names(overlap_ls) = params

for (currentParam in params)
	overlap_ls[[currentParam]] = overlap_fct(polygons_ls[[currentParam]], n_runs)

overlap_dt = rbindlist(overlap_ls, idcol = "parameter")

saveRDS(overlap_dt, paste0(path, "overlap_dt.rds"))
