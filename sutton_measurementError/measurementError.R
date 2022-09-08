
#### Aim of prog: Investigate on the measurement error
## Explanations
# The data used in this script are from the lab Integrative Ecology (PI: Dominique GRAVEL, https://ielab.recherche.usherbrooke.ca)
#
# 945 trees have been measured by two people in Sutton Canada. It was done this way:
#	1. The first person measure many trees (more than a hundred in a row) while the second person takes notes
#	2. Then, the roles are reversed: person 2 measures and person 1 takes notes.
# There should have enough measurements in a row to make sure that the second person is not influenced by the first
# This was done the last week of May 2016
# I also got data from previous campaign in Sutton from 2011 to 2013, and from 2011 to 2016
# 
## Names column data:
#	- Arbre = tree id number
#	- Esp = species
#	- Etat = state (V = alive, there are only living trees)
#	- Multi = is it a multi-trunk or not. N for No
#	- The last two columns are the values collected by the two measurers. Note that the names of the measurers appeared in the data provided
#		by the lab Gravel (Sherbrooke, Canada). To preserve the anonymity, these names were replaced by dbh1_in_mm, dbh2_in_mm
#	- The two persons measured the diameters of trees in mm at the height 1.37m (breast height).
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(MetBrewer)
library(stringi)

#### Tool functions
## Related to checking/cleaning data
# To check that each tree has only one species associated (i.e., individuals do not change species)
checkSpecies = function(dt, usedKey = "Arbre", option = "simplify")
{
	if ("usedKey" %in% names(dt))
		warning("Data table might be confused by the argument usedKey. It is both an argument of the function and a column name")

	dt[, species_nb := length(unique(Esp)), by = usedKey]
	if (dt[species_nb > 1, .N] != 0)
	{
		warning("Some trees have more than one species!")
		if (option == "simplify")
		{
			warning("You chose the option to simplify, species will be rewritten as unknown")
			dt[species_nb > 1, Esp := "unknown"]
		}
	}
	dt[, species_nb := NULL]
}

# To check that each tree id is unique (note that multi-trunk trees can have more than one row. Use the ... argument with toKeep)
checkIndividuals = function(dt, usedKey = "Arbre", option = "remove", ...)
{
	if ("usedKey" %in% names(dt))
		warning("Data table might be confused by the argument usedKey. It is both an argument of the function and a column name")
	
	dt[, nb_indiv := .N, by = usedKey]
	if (dt[, max(nb_indiv) != 1])
	{
		warning("Some individuals have more than one row of data. First simplification with function `unique'")
		dt = unique(dt)
		dt[, nb_indiv := .N, by = usedKey]
		
		if (dt[, max(nb_indiv) != 1])
		{
			warning("Despite applying unique, there are still some individuals with at least two rows of data")
			if (option == "remove")
			{
				warning("You chose the option to remove the concerned rows (except those that are to kept, set-up by the key)")
				dt = dt[nb_indiv == 1]
			}
		}
	}
	dt[, nb_indiv := NULL]
	return (dt) # No other choice than returning because I am removing some data. Not possible to do it by reference
}

# To plot traces of mcmcs from cmdstan fit objects
lazyTrace = function(draws, ...)
{
	if (!all(class(draws) %in% c("draws_array", "draws", "array")))
		stop("The class of draws should be compatible with stan (draws_array, draws, array)")
	
	n_chains = dim(draws)[2]
	n_iter = dim(draws)[1]
	cols = c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")

	min_val = min(draws)
	max_val = max(draws)
	plot(0, pch = "", xlim = c(0, n_iter), ylim = c(min_val, max_val), axes = TRUE,
		xlab = "", ylab = "", bg = "transparent")

	for (chain in 1:n_chains)
		lines(1:n_iter, draws[, chain,], type = "l", col = cols[chain])	
	
	providedArgs = list(...)
	nbArgs = length(providedArgs)

	if (nbArgs != 0)
	{
		ls_names = names(providedArgs)
		val_ind = stri_detect(str = ls_names, regex = "val")

		if (any(val_ind))
		{
			for (val in ls_names[val_ind])
				abline(h = providedArgs[[val]], col = "#CD212A")
		}
	}
}

#### Loading and cleaning data
## Load trees measured twice the same day in May 2016
treeData_remeasured = readRDS("./trees_remeasured.rds")
print(paste("The data set contains", treeData_remeasured[, .N], "measures"))
print(paste("The data set contains", length(treeData_remeasured[, unique(Arbre)]), "individuals"))

# Checking species
checkSpecies(treeData_remeasured)

# Checking uniqueness of the individuals (see crash test zone to keep few more data). I will remove the multi-trunk trees
multiTrunk = treeData_remeasured[Multi == "O", Arbre] # List the trees that have multi-trunk! Those trees for sure have more than one line

treeData_remeasured = treeData_remeasured[!(Arbre %in% multiTrunk)]

treeData_remeasured = checkIndividuals(treeData_remeasured, usedKey = "Arbre")

## Keep only living trees
treeData_remeasured = treeData_remeasured[Etat == "V"]

## Last check-up
if (treeData_remeasured[, length(unique(Arbre))] != treeData_remeasured[, .N])
	stop("Check the uniqueness of the data again!")

## Keep only the column of interest
treeData = treeData_remeasured[, .(Arbre, Esp, dbh1_in_mm, dbh2_in_mm)]

treeData[, lapply(.SD, mean), .SDcols = c("dbh1_in_mm", "dbh2_in_mm")]
treeData[, lapply(.SD, sd), .SDcols = c("dbh1_in_mm", "dbh2_in_mm")]
treeData[, lapply(.SD, range), .SDcols = c("dbh1_in_mm", "dbh2_in_mm")]

plot(treeData[, dbh1_in_mm], treeData[, dbh2_in_mm], pch = 20)

## Change colnames
setnames(treeData, c("Arbre", "Esp"), c("tree_id", "species"))

#### Estimate the routine and extreme observation errors
## Common variables
set.seed(1969-08-18) # Woodstock seed
maxIter = 3000
n_chains = 3

## Data
stanData = list(
	# Number of data
	n_trees = treeData[, .N], # Number of trees

	# Data
	dbh1 = treeData[, dbh1_in_mm],
	dbh2 = treeData[, dbh1_in_mm]
)

## Compile model
model = cmdstan_model("./measurementError.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE,
	max_treedepth = 14, adapt_delta = 0.8)

## Check-up
results$cmdstan_diagnose()
results$print()

## Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = "./", basename = paste0("estimationError-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0("./", "estimationError-", time_ended, ".rds"))
