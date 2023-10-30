
#### Aim of prog: Plot posterior extents for quadratic terms

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)

#### Source tool functions
source("./toolFunctions.R")

#### Runs
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
run = 1

ssm_list = vector(mode = "list", length = length(ls_species))
classic_list = vector(mode = "list", length = length(ls_species))
names(ssm_list) = ls_species
names(classic_list) = ls_species

## Results
for (species in ls_species)
{
	tree_path = paste0("./", species, "/")
	info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
	ssm_list[[species]] = readRDS(paste0(tree_path, info_lastRun[["file"]]))

	info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
	classic_list[[species]] = readRDS(paste0(tree_path, info_lastRun[["file"]]))

	print(paste(species, "loaded"))
}

## Plots
extentPosterior(ssm_list = ssm_list, classic_list = classic_list, ext = ".tex")
