
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)


#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Create the cluster #! USE A BASH SCRIPT TO RUN THIS R SCRIPT, EXPORT THE ARRAY_ID
# array_id = 46 #! REMOVE THIS LINE WHEN DONE WITH TEST. FAGUS_SYLVATICA
array_id = 71 #! REMOVE THIS LINE WHEN DONE WITH TEST. PICEA_ABIES
# array_id = 85 #! REMOVE THIS LINE WHEN DONE WITH TEST. PINUS_SYLVESTRIS
# array_id = 150 #! REMOVE THIS LINE WHEN DONE WITH TEST. TILIA PLATYPHYLLOS
# array_id = as.integer(Sys.getenv("ARRAY_ID")) # Each array takes care of one species
# print(paste0("array id = ", array_id))

#### Tool function
## Initiate Y_gen with reasonable values (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("results", "stanFile", "stanData")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide results, stanFile, and stanData")

	results = providedArgs[["results"]]
	stanFile = providedArgs[["stanFile"]]
	stanData = providedArgs[["stanData"]]

	n_chains = results$num_chains()

	if (!(type %in% c("mean", "median")))
		stop ("Unknown type. So far, only mean or median")

	model_init = cmdstan_model(stanFile)

	initValues = model_init$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

	list_init = vector(mode = "list", length = n_chains)
	
	index = sample(x = 1:results$metadata()$iter_sampling, size = n_chains, replace = FALSE)
	print("Formatting generated data")

	latent_growth_init_array = initValues$draws(variables = "latent_growth_init") # iter x n_chains x n_latent_growth

	for (i in 1:n_chains)
	{
		print(paste("Starting iteration", i))
		list_init[[i]] = list(
			latent_dbh_parents = as.vector(initValues$draws("latent_dbh_parents_init")[index[i], i, ]),
			latent_growth = as.vector(latent_growth_init_array[index[i], i, ]),

			potentialGrowth = as.numeric(results$draws("potentialGrowth")[index[i], i, ]),
			dbh_slope = as.numeric(results$draws("dbh_slope")[index[i], i, ]),

			pr_slope = as.numeric(results$draws("pr_slope")[index[i], i, ]),
			pr_slope2 = as.numeric(results$draws("pr_slope2")[index[i], i, ]),
			tas_slope = as.numeric(results$draws("tas_slope")[index[i], i, ]),
			tas_slope2 = as.numeric(results$draws("tas_slope2")[index[i], i, ]),

			competition_slope = as.numeric(results$draws("competition_slope")[index[i], i, ]),

			processError = as.numeric(results$draws("processError")[index[i], i, ]),
			measureError = as.numeric(results$draws("measureError")[index[i], i, ]))

		print(paste0(round(i*100/n_chains), "% done"))
	}

	return(list_init)
}

## Get name of the last run
getLastRun = function(path, begin = "growth-", extension = ".rds", format = "ymd", run = 1, getAll = FALSE)
{
	if (format != "ymd")
		stop("Format is not recognised. For now only year-month-day alias ymd works")
	
	if (!is.null(run))
		begin = paste0(begin, "run=", run, "-")
	
	ls_files = list.files(path = path, pattern = paste0("^", begin, ".*", extension, "$"))
	ls_files_split = stri_split(stri_sub(ls_files, from = stri_locate(ls_files, regex = begin)[, "end"] + 1),
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
## Paths
mainFolder = "/home/amael/project_ssm/inventories/FR IFN/processed data/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

clim_folder = "/home/amael/project_ssm/climateData/Chelsa/yearlyAverage/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((array_id < 1) | (array_id > length(ls_species)))
	stop(paste0("Array id = ", array_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

species = ls_species[array_id]
print(paste("Script running for species:", species))
treeData = treeData[speciesName_sci == species]
savingPath = paste0("./", species, "/")
if (!dir.exists(savingPath))
	stop(paste0("The species `", species, "' has no pre-run"))

## Climate
climate = readRDS(paste0(clim_folder, "FR_reshaped_climate.rds"))

## Set-up indices
indices = readRDS(paste0(mainFolder, species, "_indices.rds"))

# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1], type := "parent"]

# Compute the number of growing years per individual
indices[, nbYearsGrowth := max(year) - min(year), by = .(pointInventory_id, tree_id)]

checkUp = all(indices[, nbYearsGrowth == index_clim_end - index_clim_start])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

if (indices[, .N] != treeData[, .N])
	stop(paste0("Dimension mismatch between indices and treeData for species `", species, "`"))

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

n_indiv = unique(treeData[, .(tree_id, pointInventory_id)])[, .N]
print(paste("Number of individuals:", n_indiv))

nbYearsGrowth = unique(indices[, .(tree_id, pointInventory_id, nbYearsGrowth)])[, nbYearsGrowth]
if (length(nbYearsGrowth) != n_indiv)
	stop("Dimension mismatch between nbYearsGrowth and n_indiv")

n_hiddenState = indices[.N, index_gen]
print(paste("Number of latent dbh:", n_hiddenState))
n_latentGrowth = n_hiddenState - n_indiv

parents_index = treeData[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(pointInventory_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(pointInventory_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

#### Compute and save normalising constantes
climate_mu_sd = readRDS(paste0(savingPath, "climate_normalisation.rds"))

#### Stan model
## Define stan variables
# Common variables
maxIter = 1500
n_chains = 3
type = "median"

# Read pre-run
info_lastRun = getLastRun(path = savingPath, run = 1)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]

results = readRDS(paste0(savingPath, lastRun))
stanData_firstRun = readRDS("Picea_abies/1_stanData.rds")

## Stan data
# Data to provide to last run
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_climate = climate[, .N], # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_latentGrowth = n_latentGrowth, # Dimension of the state space vector for latent growth
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsGrowth = nbYearsGrowth, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	children_index = children_index, # Index of children in the observed data
	climate_index = indices[type == "parent", index_clim_start], # Index of the climate associated to each parent

	# Observations
	Yobs = treeData[, dbh],
	sd_dbh = treeData[, sd(dbh)],

	# Explanatory variables
	precip = climate[, pr], # Annual precipitations (sum over 12 months)
	pr_mu = climate_mu_sd[variable == "pr", mu],
	pr_sd = climate_mu_sd[variable == "pr", sd],

	tas = climate[, tas], # Annual average temperature (average over 12 months)
	tas_mu = climate_mu_sd[variable == "tas", mu],
	tas_sd = climate_mu_sd[variable == "tas", sd],

	totalTreeWeight = treeData[parents_index, totalTreeWeight] # Includes also the weight of other species!
)

# Data to provide to initiate the states
stanData_generate = stanData

stanData_generate$n_indiv_new = n_indiv
stanData_generate$n_indiv = stanData_firstRun$n_indiv

stanData_generate$n_latentGrowth_new = n_hiddenState - n_indiv
stanData_generate$n_latentGrowth = stanData_firstRun$n_latentGrowth

# Add prior informations to 
stanData$potentialGrowthPrior_mean = ifelse(type == "median",
	median(results$draws("potentialGrowth")),
	mean(results$draws("potentialGrowth")))
stanData$potentialGrowthPrior_sd = sd(results$draws("potentialGrowth"))


stanData$dbh_slopePrior_mean = ifelse(type == "median",
	median(results$draws("dbh_slope")),
	mean(results$draws("dbh_slope")))
stanData$dbh_slopePrior_sd = sd(results$draws("dbh_slope"))


stanData$pr_slopePrior_mean = ifelse(type == "median",
	median(results$draws("pr_slope")),
	mean(results$draws("pr_slope")))
stanData$pr_slopePrior_sd = sd(results$draws("pr_slope"))


stanData$pr_slope2Prior_mean = ifelse(type == "median",
	median(results$draws("pr_slope2")),
	mean(results$draws("pr_slope2")))
stanData$pr_slope2Prior_sd = sd(results$draws("pr_slope2"))


stanData$tas_slopePrior_mean = ifelse(type == "median",
	median(results$draws("tas_slope")),
	mean(results$draws("tas_slope")))
stanData$tas_slopePrior_sd = sd(results$draws("tas_slope"))


stanData$tas_slope2Prior_mean = ifelse(type == "median",
	median(results$draws("tas_slope2")),
	mean(results$draws("tas_slope2")))
stanData$tas_slope2Prior_sd = sd(results$draws("tas_slope2"))


stanData$competition_slopePrior_mean = ifelse(type == "median",
	median(results$draws("competition_slope")),
	mean(results$draws("competition_slope")))
stanData$competition_slopePrior_sd = sd(results$draws("competition_slope"))


stanData$processErrorPrior_mean = ifelse(type == "median",
	median(results$draws("processError")),
	mean(results$draws("processError")))
stanData$processErrorPrior_sd = sd(results$draws("processError"))


stanData$measureErrorPrior_mean = ifelse(type == "median",
	median(results$draws("measureError")),
	mean(results$draws("measureError")))
stanData$measureErrorPrior_sd = sd(results$draws("measureError"))

## Create initial values
initValues = init_fun(results = results, stanFile = "generate_initVals.stan", stanData = stanData_generate)

## Compile model
model = cmdstan_model("./restartGrowth.stan")

## Run model
results2 = model$sample(data = stanData, parallel_chains = n_chains, refresh = 5, chains = n_chains,
	iter_warmup = round(maxIter/3), iter_sampling = round(2*maxIter/3), adapt_delta = 0.9, max_treedepth = 13,
	# adapt_engaged = TRUE,
	step_size = results$metadata()$step_size_adaptation,
	init = initValues)


time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results2$save_output_files(dir = savingPath, basename = paste0("growth-", time_ended), timestamp = FALSE, random = TRUE)
results2$save_object(file = paste0(savingPath, "growth-", time_ended, ".rds"))

results2$cmdstan_diagnose()

results2$print(c("potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"competition_slope", "measureError", "processError"))
