
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
		
		for (j in seq_len(results_dt[, .N])) # Loop among the combinations with 'nbElements' runs
		{
			selected_runs = paste0("run_", combinations_ls[, j])
			currentCombination = paste(combinations_ls[, j], collapse = "_")

			union_polygon = polygons_ls[[selected_runs[1]]]
			inter_polygon = polygons_ls[[selected_runs[1]]]

			for (currentRun in selected_runs[2:length(selected_runs)]) # Loop among the selected runs
			{
				union_polygon = terra::union(union_polygon, polygons_ls[[currentRun]])
				union_polygon = terra::aggregate(union_polygon)
				inter_polygon = terra::intersect(inter_polygon, polygons_ls[[currentRun]])
			}
			results_dt[j, combination := currentCombination]
			area = suppressWarnings(ifelse(is.polygons(inter_polygon), terra::expanse(inter_polygon), 0))
			suppressWarnings(results_dt[j, overlap := area/terra::expanse(union_polygon)])
		}
		overlap_ls[[paste0(nbElements, "_elements")]] = results_dt
	}

	overlap_dt = rbindlist(overlap_ls, idcol = "elements")
	overlap_dt[, elements := stri_sub(str = elements, to = stri_locate_first(str = elements, regex = "_")[, "start"] - 1)]
	
	return(overlap_dt)
}

## Function to plot growth posteriors for all runs and compute their overlaps
lazyPosteriorGrowth = function(draws, printPlot = FALSE)
{
	n_runs = length(draws)

	# Get posterior
	density_from_draws = vector(mode = "list", length = n_runs)
	x = vector(mode = "list", length = n_runs)
	y = vector(mode = "list", length = n_runs)
	polygons_ls = vector(mode = "list", length = n_runs)
	names(polygons_ls) = paste0("run_", 1:n_runs)

	for (run in seq_len(n_runs))
	{
		density_from_draws[[run]] = density(sd_dbh*draws[[run]], n = 512)
		x[[run]] = density_from_draws[[run]]$x
		y[[run]] = density_from_draws[[run]]$y

		temporary = cbind(id = 1, part = 1, x[[run]], y[[run]])
		polygons_ls[[run]] = vect(temporary, type = "polygons")
	}
	min_x = min(sapply(x, min))
	max_x = max(sapply(x, max))
	max_y = max(sapply(y, max))

	min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
	max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
	
	# Plot posterior and compute overlap
	if (printPlot)
	{
		colours = MetBrewer::met.brewer("Hokusai3", n_runs)
		colours_str = grDevices::colorRampPalette(colours)(n_runs)
		colours_str_pol = paste0(colours_str, "66")
		plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", main = "", xlab = "")
		for (i in seq_len(n_runs))
		{
			lines(x = density_from_draws[[i]]$x, y = density_from_draws[[i]]$y, col = colours_str[i], lwd = 2)
			polygon(density_from_draws[[i]], col = colours_str_pol[i])
		}
	}

	overlap = overlap_fct(polygons_ls = polygons_ls, n_runs = n_runs)[elements == n_runs, overlap]
	return (overlap)
}

## Other functions
source("./toolFunctions.R")

#### Common variables
n_runs = 4 # Number of runs used in growth_subsample.R

## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

threshold_indiv = 12000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2023/01/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_de-fr-sw_12000_main.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

species = infoSpecies[1, speciesName_sci]
for (species in infoSpecies[, speciesName_sci])
{
	path = paste0("./", species, "/")

	if (infoSpecies[species, n_indiv] < threshold_indiv)
		next

	if (file.exists(paste0(path, "overlap_states.rds")))
	{
		print(paste0("Species <", species, "> already done"))
		next
	}

	results_ls = vector(mode = "list", length = n_runs)
	draws_ls = vector(mode = "list", length = n_runs)
	hiddenStates_ls = vector(mode = "list", length = n_runs)

	treeData_ls = vector(mode = "list", length = n_runs)

	sd_dbh = readRDS(paste0(path, "1_stanData.rds"))$sd_dbh # Same regardless of the run (see growth_subsample)
	indices_ls = vector(mode = "list", length = n_runs)

	oldw = getOption("warn")
	options(warn = -1)
	for (run in seq_len(n_runs))
	{
		info_lastRun = getLastRun(path = path, extension = "_main.rds$", run = run)
		lastRun = info_lastRun[["file"]]
		results_ls[[run]] = readRDS(paste0(path, lastRun))
		draws_ls[[run]] = results_ls[[run]]$draws("latent_growth")

		treeData_ls[[run]] = readRDS(paste0(path, run, "_treeData.rds"))
		indices_ls[[run]] = readRDS(paste0(path, run, "_indices.rds"))[["indices"]][, .(plot_id, tree_id, year,
			index_latent_growth, nbYearsGrowth)]
	}
	options(warn = oldw)

	treeData_dt = rbindlist(treeData_ls, idcol = "run")
	indices_dt = rbindlist(indices_ls, idcol = "run")

	treeData_dt[, occ_indiv := .N, by = .(plot_id, tree_id, year)]
	treeData = unique(treeData_dt[occ_indiv == n_runs, .(plot_id, tree_id, year, run)])

	test = data.table::merge.data.table(x = treeData, y = indices_dt, by = c("plot_id", "tree_id", "year", "run"))
	setkey(test, plot_id, tree_id, run, year) # Year is not that useful, but ensure that it is ordered

	if (test[, .N] == 0)
	{
		warning(paste0("No common individuals within ", n_runs, " runs for species <", species, ">"))
		next
	}

	setkey(test, plot_id, tree_id, run, year) # Year is not that useful, but ensure that it is ordered

	ls_plots = unique(test[, plot_id])

	overlap = numeric(length = sum(unique(test[, .(plot_id, tree_id, nbYearsGrowth)])[, nbYearsGrowth]))
	n_overlap = length(overlap)
	count = 1

	for (plot in ls_plots)
	{
		for (tree in unique(test[.(plot), tree_id]))
		{
			nb_states = unique(test[.(plot, tree), nbYearsGrowth])
			if (length(nb_states) != 1)
				stop("Non-unique number of state for a given combination plot/tree")

			for (currentState in seq_len(nb_states))
			{
				for (rr in seq_len(n_runs))
				{
					index = test[.(plot, tree, rr), ][1, index_latent_growth] + currentState - 1

					hiddenStates_ls[[rr]] = draws_ls[[rr]][, , index]
				}
				overlap[count] = lazyPosteriorGrowth(draws = hiddenStates_ls, printPlot = FALSE)

				count = count + 1
			}
		}
		print(paste0(round((count - 1)*100/n_overlap, 2), "% done"))
	}

	saveRDS(overlap, paste0(path, "overlap_states.rds"))
}
