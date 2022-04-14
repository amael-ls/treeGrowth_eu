
#### Aim of script: Prepare the data for growth (indices for multiple countries)
## Comments
# There are two parts in this program:
#	1. Standardisation of the growth data among country
#	2. Extract climate data and check that all the rasters have the same grid:
#		- It is important that all the rasters share the same grid because I am using indices in stan and I need to make sure that indices
#			and geographic locations are equivalent (i.e., an index designates the same location regardless of the raster)
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#? --------------------------------------------------------------------------------------------------------
######## PART I: Standardisation of the growth data among country
#? --------------------------------------------------------------------------------------------------------
#### Create tables for standardisations
## Table for standardising column names and names per NFI
std_colnames = data.table(standardised = c("plot_id", "tree_id", "speciesName_sci", "year", "dbh", "x", "y", "standBasalArea"),
	france = c("pointInventory_id", "tree_id", "speciesName_sci", "year", "dbh", "xLambert93", "yLambert93", "standBasalArea"),
	germany = c("plot_id", "tree_id", "speciesName_sci", "date", "dbh", "x_epsg31467", "y_epsg31467", "standBasalArea"))

## Table for standardising projection systems
std_epsg = data.table(standardised = "EPSG:4326",
	france = "EPSG:2154",
	germany = "EPSG:31467")

#### Load and reshape data
## List processed data
inputDataFolder = "/bigdata/Inventories/"
outputFolder = "/home/amael/project_ssm/inventories/"

if (!dir.exists(inputDataFolder))
	stop("Input data folder does not exist")

if (!dir.exists(outputFolder))
	dir.create(outputFolder)

NFI_ls = data.table(NFI = c("AT OWI", "CH LFI", "DE BaySF", "DE BWI", "ES IFN", "EU ICP", "FI NFI", "FR IFN", "SE NFI", "SK NIML"),
	country = c("", "", "germany_bavaria", "germany", "spain", "", "", "france", "sweden", ""))

NFI_ls[, isProcessed := "processed data" %in% list.dirs(paste0(inputDataFolder, NFI), full.name = FALSE, recursive = FALSE), by = NFI]

## Load processed data, keep only columns of interest, project to epsg 4326 (wgs84)
nb_processed = sum(NFI_ls[, isProcessed])

data_ls = vector(mode = "list", length = nb_processed)
names(data_ls) = NFI_ls[(isProcessed), NFI]

for (i in 1:nb_processed)
{
	current_country = NFI_ls[(isProcessed)][i, country]
	dataNFI = NFI_ls[(isProcessed)][i, NFI]
	data_ls[[dataNFI]] = readRDS(paste0(inputDataFolder, dataNFI, "/processed data/growth_data.rds"))

	# Subset and keep only cols
	cols = unname(unlist(std_colnames[, ..current_country]))
	data_ls[[dataNFI]] = data_ls[[dataNFI]][, ..cols]

	# Standardised names
	setnames(data_ls[[dataNFI]], old = cols, new = std_colnames[, standardised])

	# Standardised projections
	coords = vect(data_ls[[dataNFI]][, .(x, y)], geom = c("x", "y"), crs = unlist(std_epsg[, ..current_country]))
	coords = project(coords, std_epsg[, standardised])
	data_ls[[dataNFI]][, c("x", "y") := crds(coords, df = TRUE)]

	# Standardised year
	if (class(data_ls[[dataNFI]][, year]) == "Date")
		data_ls[[dataNFI]][, year := as.integer(format(year, "%Y"))]

	print(paste(current_country, "done"))
}

#### Bind data, add country to plot_id, and save data
## Binding
growthData = rbindlist(data_ls, idcol = "nfi_id", use.names = TRUE)

## Add country (note that I remain neutral regarding the status of Bavaria in Germany, although I called the column country...)
growthData[, country := NFI_ls[NFI == nfi_id, country], by = nfi_id]
growthData[, plot_id := paste(plot_id, country, sep = "_")]

## Save data
saveRDS(growthData, paste0(outputFolder, "standardised_european_growth_data.rds"))



#? --------------------------------------------------------------------------------------------------------
######## PART II: Extract climate data and check that all the rasters have the same grid
#? --------------------------------------------------------------------------------------------------------
#### Tool functions
## Function to check that the selected climate variables all cover the same timespan
timeCover_climate = function(clim_folder, selected_climate = dir(clim_folder))
{
	available_climate = dir(clim_folder)
	available_climate = available_climate[available_climate %in% selected_climate]
	if (length(available_climate) == 0)
		stop(paste0("No climate variables found among: \n- ", paste0(selected_climate, collapse = "\n- ")))

	ls_rasters_ref = sort(dir(paste0(clim_folder, available_climate[1], "/"), pattern = ".tif$"))
	checked_climate = character(length(available_climate))
	count = 1

	for (clim_var in available_climate)
	{
		ls_rasters = sort(dir(paste0(clim_folder, clim_var, "/"), pattern = ".tif$"))
		available_climate_years = as.integer(stri_sub(ls_rasters, to = stri_locate(ls_rasters, regex = ".tif")[, "start"] - 1))
		checked_climate[count] = clim_var
		
		if (!isTRUE(all.equal(ls_rasters, ls_rasters_ref)))
			return (list(sameTimeCover = FALSE, available_climate_years = available_climate_years, checked_climate = checked_climate))

		count = count + 1
	}

	return (list(sameTimeCover = TRUE, available_climate_years = available_climate_years, checked_climate = checked_climate))
}

## Function to check rasters' projection, origin, and extent
checkRasters = function(selectedVariables, years)
{
	nb_clim_variables = length(selectedVariables)
	ref_raster_ls = vector(mode = "list", length = nb_clim_variables)

	mismatch_within = data.table(variable = character((length(years) - 1)*nb_clim_variables),
		file = character((length(years) - 1)*nb_clim_variables))

	count_within_mismatch = 0
	iter = 0

	same_proj_within_var = FALSE
	same_proj_among_var = FALSE
	is_prevalent_proj = FALSE

	for (clim_var in selectedVariables)
	{
		iter = iter + 1
		ref_raster_ls[[iter]] = rast(paste0(clim_folder, clim_var, "/", years[1], ".tif"))
		ref_proj = crs(ref_raster_ls[[iter]], describe = TRUE, proj = TRUE)
		ref_origin = origin(ref_raster_ls[[iter]])
		ref_extent = ext(ref_raster_ls[[iter]])

		for (year in years[2:length(years)])
		{
			clim_rs = rast(paste0(clim_folder, clim_var, "/", year, ".tif"))
			proj = crs(clim_rs, describe = TRUE, proj = TRUE)
			check_proj = isTRUE(all.equal(ref_proj, proj))

			origin = origin(clim_rs)
			check_origin = isTRUE(all.equal(ref_origin, origin))

			extent = ext(clim_rs)
			check_extent = isTRUE(all.equal(ref_extent, extent))

			if (!check_proj | !check_origin | !check_extent)
			{
				count_within_mismatch = count_within_mismatch + 1
				mismatch_within[count_within_mismatch, c("variable", "file") := .(clim_var, year)]
			}
		}
	}

	mismatch_within = na.omit(mismatch_within)
	mismatch_within = mismatch_within[variable != ""]

	if (count_within_mismatch == 0)
		same_proj_within_var = TRUE
	
	ls_crs = vector(mode = "list", length = nb_clim_variables)
	for (iter in 1:nb_clim_variables)
		ls_crs[[iter]] = crs(ref_raster_ls[[iter]], describe = TRUE, proj = TRUE)
	ls_crs = rbindlist(ls_crs)

	ls_crs[, occurence := .N, by = proj]
	ls_crs[, variable := selectedVariables]

	if (length(ls_crs[, unique(proj)]) == 1) # To indicate if there is one most prevalent proj when projs differs
		same_proj_among_var = TRUE

	most_common_proj = unique(ls_crs[occurence == max(occurence), proj])

	if (length(most_common_proj) == 1)
		is_prevalent_proj = TRUE
	
	if (!same_proj_within_var)
		return (list(same_proj_within_var = same_proj_within_var, same_proj_among_var = same_proj_among_var,
			most_common_proj = most_common_proj, is_prevalent_proj = is_prevalent_proj, mismatch_within = mismatch_within, ls_crs = ls_crs))

	if (!same_proj_among_var)
		return (list(same_proj_within_var = same_proj_within_var, same_proj_among_var = same_proj_among_var,
			most_common_proj = most_common_proj, is_prevalent_proj = is_prevalent_proj, ls_crs = ls_crs))

	return (list(same_proj_within_var = same_proj_within_var, same_proj_among_var = same_proj_among_var,
		most_common_proj = most_common_proj, is_prevalent_proj = is_prevalent_proj))
}

## Process the mismatch data table returned by checkRaster
processMismatch = function(mismatch)
{
	should_raster_be_checked = FALSE
	ls_names = c("same_proj_within_var", "same_proj_among_var", "most_common_proj", "is_prevalent_proj")
	if (!all(ls_names %in% names(mismatch)))
		stop("From processMismatch: some data are missing in mismatch list")
	
	if (!mismatch[["same_proj_within_var"]])
		stop("Projections within variables mismatch. Check the returned data table `mismatch_within`. CRS details are provided in `ls_crs`")

	if (!mismatch[["same_proj_among_var"]])
	{
		if (!("ls_crs" %in% names(mismatch)))
			stop("From processMismatch: `ls_crs`is missing in mismatch list")

		writeLines("Projections among variables mismatch!")
		if (!mismatch[["is_prevalent_proj"]])
			stop("There is no prevalent projection. CRS details are provided in `ls_crs`")
		
		if (mismatch[["is_prevalent_proj"]])
		{
			should_raster_be_checked = TRUE
			ls_reproj_var = mismatch[["ls_crs"]][occurence != max(occurence), variable]
			writeLines(paste0("There is a prevalent projection. Reprojecting the following variables into the prevalent projection:\n- ",
				paste0(ls_reproj_var, collapse = "\n- ")))

			ref_variable = mismatch[["ls_crs"]][occurence == max(occurence), variable][1]
			ref_raster = rast(paste0(clim_folder, ref_variable, "/", timespan[["available_climate_years"]][1], ".tif"))
			
			for (clim_var in ls_reproj_var)
			{
				for (year in timespan[["available_climate_years"]])
				{
					raster_file = paste0(clim_folder, clim_var, "/", year, ".tif")
					temp = rast(raster_file)
					temp = project(x = temp, y = ref_raster)
					writeRaster(x = temp, filename = raster_file, overwrite = TRUE)
				}
			}
		}
	}
	return (should_raster_be_checked)
}

## Get all the years within a time range
getYears = function(time_range)
{
	years = stri_split(str = time_range, regex = "-", simplify = TRUE)
	year_1 = as.integer(years[1])
	year_2 = as.integer(years[2])

	return (as.character((year_1:year_2)))
}

#### Load data
## Paths
treeFolder = "/home/amael/project_ssm/inventories/"
clim_folder = "/bigdata/Predictors/climateChelsa/"

## Common variables
# Climate
selectedVariables = c("pr", "tas", "tasmin", "tasmax")
availableVariables = list.dirs(path = clim_folder, full.names = FALSE, recursive = FALSE)

if (any(!(selectedVariables %in% availableVariables)))
{
	warning(paste0("The following variables are not available! Subsetting accordingly: \n- ",
		paste0(selectedVariables[!(selectedVariables %in% availableVariables)], collapse = "\n- ")))
	
	selectedVariables = selectedVariables[selectedVariables %in% availableVariables]
}

# Tree data
treeData = readRDS(paste0(treeFolder, "standardised_european_growth_data.rds"))

timespan = timeCover_climate(clim_folder, selected_climate = c("pr", "tas", "tasmin", "tasmax"))
if (!timespan[["sameTimeCover"]])
	warning("The climate data from\n\t", clim_folder, "\ndo not all cover the same years for the tested variable. This might involve some NAs and bugs")

#### Subset tree data to cover the same timespan as climate
## Remove data from inventories that do not have associated climate data
if (!all(treeData[, unique(year)] %in% timespan[["available_climate_years"]]))
	warning(paste0("Some years of climate data are missing. The following years are removed from the tree data: \n- ",
		paste0(treeData[, unique(year)][!(treeData[, unique(year)] %in% timespan[["available_climate_years"]])], collapse = "\n- ")))

treeData = treeData[year %in% timespan[["available_climate_years"]]]
treeData[, nb := .N, by = .(nfi_id, plot_id, tree_id)]
treeData = treeData[nb > 1]
treeData[, nb := NULL]

setorder(treeData, plot_id, tree_id, year)
saveRDS(treeData, paste0(treeFolder, "standardised_european_growth_data_reshaped.rds"))

## Save time-space (i.e., time interval covered per plot)
treeCoords = unique(treeData[, .(plot_id, x, y)])
index_min = unique(treeData[treeData[, .I[year == min(year)], by = c("plot_id", "x", "y")]$V1,
	.(plot_id, x, y, year)])
index_max = unique(treeData[treeData[, .I[year == max(year)], by = c("plot_id", "x", "y")]$V1,
	.(plot_id, x, y, year)])

treeCoords[index_min, on = "plot_id", min_year := i.year]
treeCoords[index_max, on = "plot_id", max_year := i.year]
treeCoords[, unique_year_id := paste(min_year, max_year, sep = "-")]

saveRDS(treeCoords, paste0(treeFolder, "time_coordinates_growth.rds"))

#### Check raster projection, origin, and extent
mismatch = checkRasters(selectedVariables = selectedVariables, years = timespan[["available_climate_years"]])
recheckRasters = processMismatch(mismatch)
if (recheckRasters)
{
	mismatch = checkRasters(selectedVariables = selectedVariables, years = timespan[["available_climate_years"]])
	recheckRasters = processMismatch(mismatch)
	if (recheckRasters)
		stop("Problems with climate rasters, check origins and extents")	
}
message("Rasters look fine")

#### Extract data
## Compute length vector (sum of the number of years required for each plot)
length_clim = treeCoords[, sum(max_year - min_year + 1)]
combination_years = unique(treeCoords[, unique_year_id])

clim_dt = data.table()

for (j in selectedVariables)
	set(clim_dt, i = NULL, j = j, value = rep(0, length_clim))

clim_dt[, plot_id := "0"]
clim_dt[, year := 0]

count_clim_var = 0

## Stack raster, note that raster layer names cannot start by a digit, so R adds 'X'
treeYears = sort(unique(treeData[, year]))
for (clim_var in selectedVariables)
{
	clim_rs = rast(x = paste0(clim_folder, clim_var, "/", treeYears, ".tif"))
	names(clim_rs) = as.character(treeYears)
	count = 0
	start_ind = 1
	
	# Extraction
	for (time_range_id in combination_years)
	{
		count = count + 1
		selectedCoords = vect(treeCoords[unique_year_id == time_range_id, .(x, y)], crs = std_epsg[, standardised],
			geom = c("x", "y")) # Order MUST be longitude, latitude
		selectedPoint_id = treeCoords[unique_year_id == time_range_id, plot_id]
		selectedLayers = getYears(time_range_id)
		
		values = setDT(extract(x = terra::subset(x = clim_rs, subset = selectedLayers), y = selectedCoords))
		values[, plot_id := selectedPoint_id]
		values = data.table::melt(data = values, id.vars = c("plot_id"), measure.vars = selectedLayers,
			variable.name = "year", value.name = "climate", variable.factor = FALSE)
		values[, year := as.integer(year)]
		end_ind = start_ind + values[, .N] - 1
		clim_dt[start_ind:end_ind, c("plot_id", "year") := values[, .(plot_id, year)]]
		clim_dt[start_ind:end_ind, c(clim_var) := values[, climate]]

		print(paste0(round(count*100/length(combination_years)), "% done for variable ", clim_var))
		start_ind = end_ind + 1
	}
	count_clim_var = count_clim_var + 1
	print(paste0("--->", round(count_clim_var*100/length(selectedVariables)), "% total done"))
}

setorder(clim_dt, plot_id, year)

saveRDS(clim_dt, paste0(outputFolder, "europe_reshaped_climate.rds"))
