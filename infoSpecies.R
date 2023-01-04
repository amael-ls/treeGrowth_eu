
#### Aims of prog: To list some useful informations about each species of the data set and to run indices_subsample without subsampling

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)

#### Load data
treeData_folder = "/home/amael/project_ssm/inventories/growth/"
treeData = readRDS(paste0(treeData_folder, "standardised_european_growth_data_reshaped.rds"))
setindex(treeData, speciesName_sci)

#### Create info data table
## Species names and number of individuals
info = unique(treeData[, .(speciesName_sci, plot_id, tree_id)])[, .(n_indiv = .N), by = speciesName_sci]

## Order and set id which corresponds to the id used in the bash program
setkey(info, speciesName_sci)
info[, species_id := 1:.N]

## Number of plots per species
info[, n_plots := treeData[speciesName_sci, length(unique(plot_id)), on = "speciesName_sci"], by = speciesName_sci]

## Number and list of NFIs per species
info[, n_nfi := treeData[speciesName_sci, length(unique(nfi_id)), on = "speciesName_sci"], by = speciesName_sci]

info[, ls_nfi := paste(treeData[speciesName_sci, unique(nfi_id), on = "speciesName_sci"], collapse = ", "), by = speciesName_sci]
info[, ls_countries := paste(treeData[speciesName_sci, unique(country), on = "speciesName_sci"], collapse = ", "), by = speciesName_sci]

#### Save informations
saveRDS(info, "./speciesInformations.rds")

#### Run for each species the function indices_subsample with run = "full" (i.e., no subsampling)
source("./indices_subsample.R")
clim_folder = treeData_folder

for (species in info[c(1:5, 17, 48), speciesName_sci]) # Do it for all species later!
	indices_subsample("full", treeData[speciesName_sci == species], paste0("./", species, "/"), treeData_folder, clim_folder)
