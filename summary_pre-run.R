
#### Aim of prog: To sum-up the pre-runs, such as mean and sd of estimated parameters among species and within species (if run > 1)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
## Get fixed values parameters
getParams = function(model_cmdstan, params_names, type = "mean")
{
	if (!(type %in% c("mean", "median")))
		stop("Unknown type. Please choose median or mean")
	
	vals = numeric(length(params_names))
	names(vals) = params_names
	for (i in 1:length(params_names))
	{
		vals[i] = ifelse(type == "mean",
			mean(model_cmdstan$draws(params_names[i])),
			median(model_cmdstan$draws(params_names[i])))
	}
	return (vals)
}

## Get name of the last run
getLastRun = function(path, begin = "^growth-", extension = ".rds$", format = "ymd", run = NULL, getAll = FALSE, hour = TRUE)
{
	if (format != "ymd")
		stop("Format is not recognised. For now only year-month-day alias ymd works")
	
	if (!is.null(run))
	{
		print(paste("Searching among runs =", run))
		begin = paste0(begin, "run=", run, "-")
	}
	
	ls_files = list.files(path = path, pattern = paste0(begin, ".*", extension))

	if (length(ls_files) == 0)
	{
		warning(paste0("No file detected in the folder '", path, "'. You were looking for '", begin, "*", extension, "'"))
		return (list(file = NA, time_ended = NA))
	}

	if (is.null(run))
		ls_files = ls_files[!stri_detect(str = ls_files, regex = paste0(begin, "run="))]

	ls_files_split = stri_split(
		str = stri_sub(str = ls_files,
			from = stri_locate(ls_files, regex = begin)[, "end"] + 1,
			to = stri_locate_last(ls_files, regex = "_[[:digit:]].*.rds")[, "start"] - 1),
		regex = "-", simplify = TRUE)
	
	if (format == "ymd") # year month day
		dt = data.table(file = ls_files, year = ls_files_split[, 1], month = ls_files_split[, 2], day = ls_files_split[, 3])

	if (hour)
	{
		dt[, c("hour", "minute") := as.list(stri_split(str = stri_sub(str = file,
				from = stri_locate_last(file, regex = "_")[, "end"] + 1,
				to = stri_locate_last(file, regex = ".rds")[, "start"] - 1),
			regex = "h", simplify = TRUE)), by = file]
	}

	dt[stri_detect(str = day, regex = extension), day := stri_sub(str = day, to = stri_locate_first(str = day, regex = "_")[,"start"] - 1)]

	setorder(dt, year, month, day, hour, minute)
	if (getAll)
		return (list(file = dt[.N, file], time_ended = paste(dt[.N, year], dt[.N, month], dt[.N, day], sep = "-"), allFiles = dt))
	return (list(file = dt[.N, file], time_ended = paste(dt[.N, year], dt[.N, month], dt[.N, day], sep = "-")))
}

isProcessed = function(path, multi, lim_time, begin = "^growth-", extension = ".rds$", format = "ymd", lower = 1, upper = 4)
{
	if (class(lim_time) != "Date")
		lim_time = as.Date(lim_time)
	run_vec = lower:upper
	n_runs = length(run_vec)
	processed = rep(FALSE, n_runs)
	if (multi)
	{
		for (i in 1:n_runs)
		{
			run = run_vec[i]
			date_run = getLastRun(path, begin = begin, extension = extension, format = format, run = i)[["time_ended"]]
			date_run = as.Date(date_run)
			if (!is.na(date_run) & (date_run > lim_time))
				processed[i] = TRUE
		}
	} else {
		date_run = getLastRun(path, begin = begin, extension = extension, format = format, run = 1)[["time_ended"]]
		date_run = as.Date(date_run)
		if (!is.na(date_run) & (date_run > lim_time))
			processed = TRUE
	}
	return (all(processed))
}

## Function to expand the basic names when there is more than one NFI
expand = function(base_names, nb_nfi, patterns = c("Obs", "proba"))
{
	if (nb_nfi < 1)
		stop("Nothing to expand")
	
	new_names = vector(mode = "list", length = length(patterns))
	old_names = base_names
	for (i in 1:length(patterns))
	{
		reg = patterns[i]
		toModify = base_names[stri_detect(base_names, regex = reg)]
		base_names = base_names[!stri_detect(base_names, regex = reg)]
		new_names[[i]] = character(length = nb_nfi*length(toModify))
		for (j in 1:length(toModify))
			new_names[[i]][((j - 1)*nb_nfi + 1):(j*nb_nfi)] = paste0(toModify[j], "[", 1:nb_nfi, "]")
	}

	combined_names = c(base_names, unlist(new_names))
	return(list(old_names = old_names, new_names = combined_names))
}

#### Load data
## Tree data
mainFolder = "/home/amael/project_ssm/inventories/growth/"
treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))

## List species, I assume folders start by an upper case letter, followed by a lower case
ls_folders = list.dirs(path = "./", full.names = FALSE, recursive = FALSE)
ls_folders = ls_folders[stri_detect(ls_folders, regex = "^[A-Z][a-z].*")]

treeData = treeData[speciesName_sci %in% ls_folders]

ls_species = unique(treeData[, .(tree_id, plot_id, speciesName_sci)])[, .N, by = speciesName_sci]
setnames(ls_species, old = c("speciesName_sci", "N"), new = c("species", "n_indiv"))
setkey(ls_species, species)

threshold = 8000
lim_time = as.Date("2022-04-23")

ls_species[, multiRun := ifelse(n_indiv > threshold, TRUE, FALSE)]
ls_species[, processed := isProcessed(species, multiRun, lim_time, lower = 1, upper = 4), by = species]

ls_species = ls_species[(isProcessed)]

base_names = c("averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc")

## Load tree data
treeFolder = "/home/amael/project_ssm/inventories/growth/"
treeData = readRDS(paste0(treeFolder, "standardised_european_growth_data_reshaped.rds"))
treeData = treeData[speciesName_sci %in% ls_species[, species]]

# Compute number of NFIs covered per species
ls_nfi = treeData[, .(nb_nfi = length(unique(nfi_id))), by = speciesName_sci]
setnames(ls_nfi, old = "speciesName_sci", new = "species")

## Merge ls_nfi and ls_species
ls_species = merge.data.table(ls_species, ls_nfi, by = "species")

## Load results
ls_params = vector(mode = "list", length = ls_species[1:2, sum(nbResults)]) #! REMOVE THE 1:2
iter = 0
for (species in ls_species[1:2, species]) #! REMOVE THE 1:2
{
	path = paste0("./", species, "/")
	params_names = base_names
	
	if (ls_species[species, nb_nfi] > 1) #! Not used, because the priors are the same for all the NFIs, so useless to distinguish them
		params_names = expand(base_names, ls_species[species, nb_nfi])[["new_names"]]

	if (ls_species[species, multiRun])
	{
		for (run in 1:ls_species[species, nbResults])
		{
			info_lastRun = getLastRun(path = path, run = run)
			lastRun = info_lastRun[["file"]]
			results = readRDS(paste0(path, lastRun))

			iter = iter + 1
			ls_params[[iter]] = getParams(results, base_names)
		}
		next;
	}

	run = 1
	info_lastRun = getLastRun(path = path, run = run)
	lastRun = info_lastRun[["file"]]
	results = readRDS(paste0(path, lastRun))

	iter = iter + 1
	ls_params[[iter]] = getParams(results, base_names)
}

ll = ls_params[[2]]
as.data.table(ll)
ll = rbindlist(ls_params)
