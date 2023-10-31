
#### Aim of prog: Plot different visualisations of growth versus predictors (diameters, environmental factors)
## Comment
# For each 'growth vs environment' plot, the dbh is set at its optimum value (i.e., highest growth), and the other environmental variables are
#	set to 0 (i.e., their average since it is scaled)
# For the 'growth vs dbh' plot, all the environmental variables are set to 0 (i.e., their average since it is scaled)
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)
library(terra)

#### Tool functions
source("toolFunctions.R")

######## Part I: growth response to environment and to time series predictions
#### Get data
## Common variables
args = commandArgs(trailingOnly = TRUE) # Example: args = c("Betula pendula", "1", "tas", "pr", "ph", "ba")
if (length(args) < 3)
	stop("Supply the species_id, run_id, and at least one variable among pr, tas, and ph!", call. = FALSE)

(species = as.character(args[1]))
(run = as.integer(args[2]))
(variable = as.character(args[3:length(args)]))

if (!all(variable %in% c("pr", "tas", "ph", "ba")))
	stop("The variables must be among pr, tas, ph, or ba")

#### Plot residuals fit
# residuals_fit(species, run, filenamePattern = "_residuals_fit.pdf")

#### Plot growth vs variables
if ((file.exists("./speciesInformations.rds")) && (file.exists("./speciesInformations_runs.rds")))
{
	info = readRDS("./speciesInformations.rds")
	info_runs = readRDS("./speciesInformations_runs.rds")
	ls_info = list(info = info, range_subset = info_runs)
} else {
	ls_info = infoSpecies()
}

ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
info_runs = info_runs[ls_species] # speciesName_sci is the key of info_runs, so no need for %in% ls_species
info = info[ls_species] # speciesName_sci is the key of info, so no need for %in% ls_species

# All the variables range are within the species range, i.e., min(xyz_025) > min(xyz) among all species, and same with max
checkup = plotGrowth(species = species, run = run, ls_info = ls_info, variables = variable, extension = "tex", caption = FALSE,
	span = list(pr = info[, c(min(pr_025), max(pr_975))], tas = info[, c(min(tas_025), max(tas_975))],
		ph = info[, c(min(ph_025), max(ph_975))], ba = info[, c(min(ba_025), max(ba_975))]),
	pr_min = info[species, pr_025], pr_max = info[species, pr_975],
	tas_min = info[species, tas_025], tas_max = info[species, tas_975],
	ph_min = info[species, ph_025], ph_max = info[species, ph_975],
	ba_min = info[species, ba_025], ba_max = info[species, ba_975])
print(checkup)

plotGrowth_dbh(species = species, run = run, ls_info = ls_info, caption = FALSE, extension = "tex",
	dbh_min = info_runs[species, dbh_025_1], dbh_max = info_runs[species, dbh_975_1],
	span = info_runs[, c(min(dbh_025_1), max(dbh_975_1))])
