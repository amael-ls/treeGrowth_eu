
#### Aim of function: Extract climate for a given species
extractClimate = function(...)
{
	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)
	nbArgs = length(providedArgs)

	isConform = FALSE

	if (all(c("climate", "indices") %in% ls_names))
	{
		print("Using climate and indices. This extract climate for the whole species, no subsampling")
		
		climate = providedArgs[["climate"]]
		indices = providedArgs[["indices"]]

		indices = unique(indices[, .(plot_id, index_clim_start, index_clim_end)])
		indices[, nbYears := index_clim_end - index_clim_start + 1]
		tot_nbYears = indices[, sum(nbYears)]

		extract_indices = integer(tot_nbYears)
		end = 0
		for (i in 1:indices[, .N])
		{
			start = end + 1
			end = start + indices[i, nbYears] - 1
			extract_indices[start:end] = indices[i, index_clim_start]:indices[i, index_clim_end]
		}

		isConform = TRUE
	}

	if (all(c("path", "run") %in% ls_names))
	{
		print(paste("Using path and run. This extract climate for subsample", run))
		
		path = providedArgs[["path"]]
		run = providedArgs[["run"]]

		indices = readRDS(paste0(path, run, "_indices.rds"))
		stanData = readRDS(paste0(path, run, "_stanData.rds"))
		
		#! TO FINISH HERE

		isConform = TRUE
	}

	if (all(c("path_ind", "path_clim", "run") %in% ls_names))
	{
		print(paste("Using path for index and climate, and run. This extract climate for subsample", run))
		
		path = providedArgs[["path"]]
		run = providedArgs[["run"]]

		climate = readRDS(paste0(path_clim, "europe_reshaped_climate.rds"))
		indices = readRDS(paste0(path_ind, run, "_indices.rds"))

		indices = unique(indices[, .(plot_id, index_clim_start, index_clim_end)])
		indices[, nbYears := index_clim_end - index_clim_start + 1]
		tot_nbYears = indices[, sum(nbYears)]

		extract_indices = integer(tot_nbYears)
		end = 0
		for (i in 1:indices[, .N])
		{
			start = end + 1
			end = start + indices[i, nbYears] - 1
			extract_indices[start:end] = indices[i, index_clim_start]:indices[i, index_clim_end]
		}
		
		isConform = TRUE
	}

	if (!isConform)
		stop("You must provide climate/indices together or path/run together")
	


}