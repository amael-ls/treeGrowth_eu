
#### Aim of function: Extract climate for a given species
extractClimate = function(...)
{
	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)

	isConform = FALSE

	climate_activated = FALSE
	indices_activated = FALSE
	ph_activated = FALSE
	ph_alone = FALSE
	
	if ("indices" %in% ls_names)
		indices_activated = TRUE

	if ("climate" %in% ls_names)
	{
		climate_activated = TRUE
		climate = providedArgs[["climate"]]
	}

	if ("ph" %in% ls_names)
	{
		ph_activated = TRUE
		ph = providedArgs[["ph"]]
	}
	
	if (indices_activated)
	{
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

		if (climate_activated & ph_activated)
			predictors = data.table::merge.data.table(climate, ph, by = "plot_id")

		if (climate_activated & !ph_activated)
			predictors = climate

		if (!climate_activated & ph_activated)
		{
			ph = data.table::merge.data.table(ph, indices[, .(plot_id, nbYears)], by = "plot_id")
			predictors = data.table(plot_id = rep(ph[, plot_id], times = ph[, nbYears]), ph = rep(ph[, ph], times = ph[, nbYears]))
			ph_alone = TRUE
		}

		if (climate_activated | ph_activated)
			isConform = TRUE
	}

	if (all(c("path", "run") %in% ls_names))
	{
		print(paste("Using path and run. This extract climate for subsample", run))
		
		path = providedArgs[["path"]]
		run = providedArgs[["run"]]

		indices = readRDS(paste0(path, run, "_indices.rds"))
		stanData = readRDS(paste0(path, run, "_stanData.rds"))
		
		precip_name = names(stanData)[stringi::stri_detect_regex(names(stanData), "pr")]
		precip_name = precip_name[!stringi::stri_detect_regex(precip_name, "_mu")]
		precip_name = precip_name[!stringi::stri_detect_regex(precip_name, "_sd")]

		temperature_name = names(stanData)[stringi::stri_detect_regex(names(stanData), "tas")]
		temperature_name = temperature_name[!stringi::stri_detect_regex(temperature_name, "_mu")]
		temperature_name = temperature_name[!stringi::stri_detect_regex(temperature_name, "_sd")]

		ph_name = names(stanData)[stringi::stri_detect_regex(names(stanData), "ph")]
		ph_name = ph_name[!stringi::stri_detect_regex(ph_name, "_mu")]
		ph_name = ph_name[!stringi::stri_detect_regex(ph_name, "_sd")]

		predictors = data.table::data.table(temp1 = stanData[[precip_name]], temp2 = stanData[[temperature_name]],
			temp3 = stanData[[ph_name]])
		data.table::setnames(predictors, new = c(precip_name, temperature_name, ph_name))

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

	if (all(c("path_ind", "path_clim", "path_ph", "run") %in% ls_names))
	{
		print(paste("Using path for index and climate, and run. This extract climate for subsample", run))
		
		path_ind = providedArgs[["path_ind"]]
		path_clim = providedArgs[["path_clim"]]
		path_ph = providedArgs[["path_ph"]]
		run = providedArgs[["run"]]

		indices = readRDS(paste0(path_ind, run, "_indices.rds"))
		climate = readRDS(paste0(path_clim, "europe_reshaped_climate.rds"))
		ph = readRDS(paste0(path_clim, "europe_reshaped_soil.rds"))

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

		predictors = data.table::merge.data.table(climate, ph, by = "plot_id")
		
		isConform = TRUE
	}

	if (!isConform)
		stop("You must provide climate/indices together or path/run together")

	return (list(extract_indices = extract_indices, extracted_predictors = if (ph_alone) predictors else predictors[extract_indices]))
}
