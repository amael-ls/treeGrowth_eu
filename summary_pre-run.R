
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
getLastRun = function(path, begin = "growth-", extension = ".rds", format = "ymd", run = NULL, getAll = FALSE)
{
	if (format != "ymd")
		stop("Format is not recognised. For now only year-month-day alias ymd works")
	
	if (!is.null(run))
	{
		print(paste("Searching among runs =", run))
		begin = paste0(begin, "run=", run, "-")
	}
	
	ls_files = list.files(path = path, pattern = paste0("^", begin, ".*", extension, "$"))

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
	n = length(ls_files)
	
	if (format == "ymd") # year month day
		dt = data.table(file = ls_files, year = ls_files_split[, 1], month = ls_files_split[, 2], day = ls_files_split[, 3])

	dt[stri_detect(str = day, regex = extension), day := stri_sub(str = day, to = stri_locate_first(str = day, regex = "_")[,"start"] - 1)]
	setorder(dt, year, month, day)
	if (getAll)
		return (list(file = dt[.N, file], time_ended = paste(dt[.N, year], dt[.N, month], dt[.N, day], sep = "-"), allFiles = dt))
	return (list(file = dt[.N, file], time_ended = paste(dt[.N, year], dt[.N, month], dt[.N, day], sep = "-")))
}

#### Load data
## List species, I assume folders start by an upper case letter, followed by a lower case
ls_folders = list.dirs(path = "./", full.names = FALSE, recursive = FALSE)
ls_folders = ls_folders[stri_detect(ls_folders, regex = "^[A-Z][a-z].*")]

ls_species = data.table(species = ls_folders, nbResults = integer(length(ls_folders)))
setkey(ls_species, species)

for (species in ls_species[, species])
	ls_species[species, nbResults := length(list.files(path = paste0("./", species), pattern = "^growth.*.rds$"))]

ls_species[, isProcessed := ifelse(nbResults > 0, TRUE, FALSE)]
ls_species[, multiRun := ifelse(nbResults > 1, TRUE, FALSE)]

ls_species = ls_species[(isProcessed)]

## Load results
for (species in ls_species[, species])
{
	if (ls_species[species, multiProcessed])
	{
		next;
	}
	results = readRDS() 
}

