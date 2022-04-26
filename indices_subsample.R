
#### Aim of prog: match climate to tree inventories
## Explanations
# I do not join the climate with tree data. Instead I create a data table called indices. This data table is key to use rstan/rcmdstan
# after! Indeed, indices contains the indices for climate:
# 	tree x in plot y and time t is related to climate[i],
# 	tree x in plot y and time t + k is related to climate[i + k],
# The data table indices also contains the matching dbh indices in the generated data:
#	tree x is measured in 2000 and in 2005. The index in the data are for instance 4 and 5.
#	However, the yearly time-step state space model will fill the missing years, so that the year 2005 will be in position 9!

indices_subsample = function(species_id, run_id, treeData, savingPath, mainFolder, climFolder)
{
	#### Tool functions
	## Function to fill the gap between two years and get the index of the provided years
	fillYears = function(years)
	{
		if (length(years) < 2)
			stop("From fillYears: Their should be at least two years to fill the gaps")

		if (is.unsorted(years))
			stop("From fillYears: years are assumed to be sorted!")

		fill_years = years[1]:years[length(years)]
		indices = which(fill_years %in% years)
		
		return (list(fill_years = fill_years, indices = indices))
	}

	## Function to create the individual indices for the state-space model (stan)
	indices_state_space = function(trees_NFI)
	{
		count = 0
		start = 0
		end = 0
		iter = 0

		if (length(trees_NFI[, unique(speciesName_sci)]) > 1)
			stop(paste0("from indices_state_space: Their should be only one species of trees. Currently there are: \n- ",
				paste0(trees_NFI[, unique(speciesName_sci)], collapse = "\n- ")))
		
		nbIndiv = unique(trees_NFI[, .(plot_id, tree_id)])[, .N]
		length_filled_years = sum(trees_NFI[, max(year) - min(year) + 1, by = .(plot_id, tree_id)][, V1])

		indices = data.table(year = integer(length_filled_years), tree_id = integer(length_filled_years),
			plot_id = character(length_filled_years), index_gen = integer(length_filled_years),
			index_clim_start = integer(length_filled_years), index_clim_end = integer(length_filled_years))

		for (plot in trees_NFI[, unique(plot_id)])
		{
			for (indiv in trees_NFI[plot_id == plot, unique(tree_id)])
			{
				years_indices = fillYears(trees_NFI[plot_id == plot & tree_id == indiv, year])
				start = end + 1
				end = end + length(years_indices[["fill_years"]])
				indices[start:end, year := years_indices[["fill_years"]]]
				indices[start:end, tree_id := indiv]
				indices[start:end, plot_id := plot]
				indices[years_indices[["indices"]] + count, index_gen := years_indices[["indices"]] + count]
				count = count + years_indices[["indices"]][length(years_indices[["indices"]])]
				iter = iter + 1
				if (iter %% 1000 == 0)
					print(paste(round(iter*100/nbIndiv, digits = 1), "% done"))
			}
		}
		print("100% done")
		return (indices)
	}

	## Function to create the individual indices to match with climate (stan)
	# This function will modify the indices data table by reference and avoiding a copy.
	indices_climate = function(indices_dt, climate, time_space)
	{
		for (plot in indices_dt[, unique(plot_id)])
		{
			min_year = time_space[plot_id == plot, min_year]
			max_year = time_space[plot_id == plot, max_year]
			tree = indices_dt[plot_id == plot, unique(tree_id)][1]
			for (tree in indices_dt[plot_id == plot, unique(tree_id)])
			{
				clim_start = climate[plot_id == plot & year == min_year, row_id]
				clim_end = climate[plot_id == plot & year == max_year, row_id]

				indices_dt[tree_id == tree & plot_id == plot,
					c("index_clim_start", "index_clim_end") := .(clim_start, clim_end)]
			}
		}
		print("100% done")
	}

	#### Load data
	## Climate
	climate = readRDS(paste0(climFolder, "europe_reshaped_climate.rds"))
	climate[, row_id := 1:.N]

	## Time space coordinates (combination of where and whem)
	time_space = readRDS(paste0(mainFolder, "time_coordinates_growth.rds"))

	#### Create indices
	## For the state space model (*.stan file)
	indices = indices_state_space(trees_NFI = treeData)

	## For the climate
	indices_climate(indices, climate, time_space)

	## Remove the zeros #! (because of a quick and dirty solution in indices_state_space, I might optimise this code later)
	indices = indices[index_gen != 0]

	## Add an index per plot
	indices[, plot_index := .GRP, by = plot_id]

	## Add an index per country
	indices[, nfi := stri_sub(plot_id, to = stri_locate_first(plot_id, regex = "_")[, "start"] - 1)]
	indices[, nfi_index := .GRP, by = nfi]

	## Saving indices for the chosen species
	saveRDS(indices, paste0(savingPath, run_id, "_indices.rds"))
	return (indices)
}
