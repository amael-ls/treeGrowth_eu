
#### Aim of script: Prepare the data for growth (indices for multiple countries)
## Comments
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(terra)

#### Create tables for standardisations
## Table for standardising column names and names per NFI
std_colnames = data.table(standardised = c("plot_id", "tree_id", "speciesName_sci", "year", "dbh", "x", "y", "standBasalArea"),
	france = c("pointInventory_id", "tree_id", "speciesName_sci", "year", "dbh", "xLambert93", "yLambert93", "standBasalArea"),
	germany = c("plot_id", "tree_id", "speciesName_sci", "date", "dbh", "x_epsg31467", "y_epsg31467", "standBasalArea"))

## Table for standardising projection systems
std_epsg = data.table(standardised = "EPSG:4326",
	france = "EPSG:2154",
	germany = "EPSG:31467")

#### Load data
## List processed data
inputDataFolder = "/bigdata/Inventories/"
outputFolder = "/home/amael/project_ssm/inventories/"

if (!dir.exists(inputDataFolder))
	stop("Input data folder does not exist")

if (!dir.exists(outputFolder))
	dir.create(outputFolder)

NFI_ls = data.table(NFI = c("AT OWI", "CH LFI", "DE BaySF", "DE BWI", "ES IFN", "EU ICP", "FI NFI", "FR IFN", "SE NFI", "SK NIML"),
	country = c("", "", "germany/bayern", "germany", "spain", "", "", "france", "sweden", ""))

NFI_ls[, isProcessed := "processed data" %in% list.dirs(paste0(inputDataFolder, NFI), full.name = FALSE, recursive = FALSE), by = NFI]

## Load processed data, keep only columns of interest, project to epsg 4326 (wgs84)
nb_processed = sum(NFI_ls[, isProcessed])

data_ls = vector(mode = "list", length = nb_processed)
names(data_ls) = NFI_ls[(isProcessed), NFI]

i = 1
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

growthData = rbindlist(data_ls, idcol = "nfi_id", use.names = TRUE)

saveRDS(growthData, paste0(outputFolder, "standardised_european_growth_data.rds"))
