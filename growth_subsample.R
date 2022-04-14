
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)


#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(terra)

#### Get parameters for run
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
  stop("Supply the species_id, run_id, and max_indiv as command line arguments!", call. = FALSE)

# species_id = 52 #! REMOVE THIS LINE WHEN DONE WITH TEST. FAGUS_SYLVATICA

species_id = 79 #! REMOVE THIS LINE WHEN DONE WITH TEST. PICEA_ABIES
run_id = 1
max_indiv = 8000
# species_id = as.integer(args[1])
# run_id = as.integer(args[2])
# max_indiv = as.integer(args[3])

#### Tool function
## Initiate Y_gen with reasonable values (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "n_latentGrowth", "average_yearlyGrowth", "nbYearsGrowth", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh_parents, n_latentGrowth, average_yearlyGrowth, nbYearsGrowth, and normalise")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	n_latentGrowth = providedArgs[["n_latentGrowth"]]
	average_yearlyGrowth = providedArgs[["average_yearlyGrowth"]]
	nbYearsGrowth = providedArgs[["nbYearsGrowth"]]
	normalise = providedArgs[["normalise"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	if (!normalise)
		Warning("This function has been coded for normalise only. The parameters potential growth, etc are scaled with sd_dbh = 135")

	if (any(average_yearlyGrowth == 0))
	{
		warning("Some average yearly growth were 0. They have been replaced by 0.5")
		average_yearlyGrowth[average_yearlyGrowth == 0] = 0.5
	}

	n_indiv = length(dbh_parents)
	Y_gen = rgamma(n_indiv, dbh_parents^2/0.5, dbh_parents/0.5) # Average = dbh_parents, variance = 0.5

	latent_growth_gen = numeric(n_latentGrowth)
	counter_growth = 0
	for (i in 1:n_indiv) # Not that this forbid trees to shrink
	{
		for (j in 1:nbYearsGrowth[i])
		{
			counter_growth = counter_growth + 1
			latent_growth_gen[counter_growth] = rgamma(n = 1, shape = average_yearlyGrowth[i]^2/0.5, rate = average_yearlyGrowth[i]/0.5)
		}
	}

	if (any(latent_growth_gen == 0))
	{
		warning("Some generated latent growth were 0. There have been replaced by 0.1 (before standardising)")
		latent_growth_gen[latent_growth_gen == 0] = 0.1
	}
	
	# All the following initialisations are based on a previous run on Tilia platyphyllos. The set sd is 4 times the estimated sd
	potentialGrowth = rnorm(n = 1, mean = -3.98, sd = 0.08)
	dbh_slope = rnorm(n = 1, mean = 0.10, sd = 0.04)
	
	pr_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)
	pr_slope2 = rnorm(n = 1, mean = 0, sd = 0.04)
	tas_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)
	tas_slope2 = rnorm(n = 1, mean = 0, sd = 0.04)
	competition_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)

	measureError = rgamma(n = 1, shape = 0.011^2/0.0015, rate = 0.011/0.0015)
	processError = rgamma(n = 1, shape = 0.0009^2/0.00012, rate = 0.0009/0.00012)

	# Normalise dbh
	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
		latent_growth_gen = latent_growth_gen/sd_dbh
	}

	return(list(latent_dbh_parents = Y_gen, latent_growth = latent_growth_gen,
		potentialGrowth = potentialGrowth, dbh_slope = dbh_slope,
		pr_slope = pr_slope, pr_slope2 = pr_slope2,
		tas_slope = tas_slope, tas_slope2 = tas_slope2,
		competition_slope = competition_slope))
}

## Function to compute the mean and sd of given variables in a data table
normalisation = function(dt, colnames = names(df), folder = "./", filename = "normalisation.rds", rm_na = TRUE, ...)
{
	if (rm_na)
		print("Warning: rm_na is activated, normalisation won't take NA into account")

	if (!is.data.table(dt))
		stop("This function is written for data table only")

	if (stri_sub(str = folder, from = stri_length(folder)) != "/")
		folder = paste0(folder, "/")

	if (!all(colnames %in% names(dt)))
	{
		warning(paste0("The following columns do not exist in the provided data and are ignored:\n- ",
			paste0(colnames[!(colnames %in% names(dt))], collapse = "\n- ")))
		colnames = colnames[colnames %in% names(dt)]
	}

	providedArgs = list(...)
	providedArgs_names = names(providedArgs)
	if ("indices" %in% providedArgs_names)
	{
		if (!any(c("col_ind", "col_ind_start", "col_ind_end") %in% providedArgs_names))
			stop("Indices provided without any column selected")
		
		ind = providedArgs[["indices"]]
		
		if ("col_ind" %in% providedArgs_names)
		{
			col_ind = providedArgs[["col_ind"]]
			if (!(col_ind %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind))

			rowsToKeep = indices[, ..col_ind]

			if (any(c("col_ind_start", "col_ind_end") %in% providedArgs_names))
				warning("col_ind_start or col_ind_end ignored")
		}

		if (("col_ind_start" %in% providedArgs_names) & !("col_ind" %in% providedArgs_names))
		{
			if (!("col_ind_end" %in% providedArgs_names))
				stop("A starting index is provided but there is no stopping index")
			
			col_ind_start = providedArgs[["col_ind_start"]]
			if (!(col_ind_start %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind_start))

			col_ind_start = providedArgs[["col_ind_start"]]

			col_start = unique(indices[[col_ind_start]])

			col_ind_end = providedArgs[["col_ind_end"]]
			if (!(col_ind_end %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind_end))
			
			col_ind_end = providedArgs[["col_ind_end"]]

			col_end = unique(indices[[col_ind_end]])

			if (length(col_start) != length(col_end))
				stop("Starting and ending indices length mismatches")

			rowsToKeep = integer(length = sum(col_end - col_start + 1))
			count = 1
			for (i in 1:length(col_start))
			{
				rowsToKeep[count:(count + col_end[i] - col_start[i])] = col_start[i]:col_end[i]
				count = count + col_end[i] - col_start[i] + 1
			}
		}
	}
	
	n = length(colnames)
	mu_sd = data.table(variable = character(n), mu = numeric(n), sd = numeric(n))

	if (!("indices" %in% providedArgs_names))
		mu_sd[, c("variable", "mu", "sd") := .(colnames, as.matrix(dt[, lapply(.SD, mean, na.rm = rm_na), .SDcols = colnames])[1,],
			as.matrix(dt[, lapply(.SD, sd, na.rm = rm_na), .SDcols = colnames])[1,])]

	if ("indices" %in% providedArgs_names)
	{

		mu_sd[, c("variable", "mu", "sd") := .(colnames, as.matrix(dt[rowsToKeep, lapply(.SD, mean, na.rm = rm_na), .SDcols = colnames])[1,],
			as.matrix(dt[rowsToKeep, lapply(.SD, sd, na.rm = rm_na), .SDcols = colnames])[1,])]
	}

	saveRDS(mu_sd, file = paste0(folder, filename))
	print(paste0("files containing coefficients saved at: ", folder, filename))
}

## Function to subsample the data either spatially (the number of individuals might not be the targeted number) or numerically
subsampling = function(dt, n_indiv_target, mode = "spatial")
{
	if (!any(c("spatial", "numeric") %in% mode))
		stop("Unknown mode, please choose spatial or numeric")

	mean_dbh_beforeSubsample = dt[, mean(dbh)]
	sd_dbh_beforeSubsample = dt[, sd(dbh)]
	quantile_beforeSubsample_25_75 = quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))
	
	if (mode == "spatial")
	{
		n_indiv_per_plot_avg = dt[, length(unique(tree_id)), by = plot_id][, mean(V1)]
		n_plots_sampling = round(n_indiv_target/n_indiv_per_plot_avg)
		coords = unique(dt[, plot_id])
		sample_plots = sample(x = coords, size = n_plots_sampling)

		setkey(dt, plot_id)
		dt = dt[sample_plots]
		setorder(dt, plot_id, tree_id, year)
	}

	if (mode == "numeric")
	{
		setkey(dt, tree_id, plot_id)
		parents_index = dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
		sampled_indices = sort(sample(x = parents_index, size = n_indiv_target, replace = FALSE))
		chosen_individuals = dt[sampled_indices, .(tree_id, plot_id)]
		dt = dt[chosen_individuals]
		setorder(dt, plot_id, tree_id, year)
	}

	diffAverage = (dt[, mean(dbh)]/mean_dbh_beforeSubsample > 1.05) | (dt[, mean(dbh)]/mean_dbh_beforeSubsample < 0.95)
	diffSD = (dt[, sd(dbh)]/sd_dbh_beforeSubsample > 1.05) | (dt[, sd(dbh)]/sd_dbh_beforeSubsample < 0.95)
	diffQuantile_25_75 = any((quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))/quantile_beforeSubsample_25_75 > 1.05) |
		(quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))/quantile_beforeSubsample_25_75 < 0.95))

	n_indiv = unique(dt[, .(tree_id, plot_id)])[, .N]

	return(list(sampledData = dt, diffAverage = diffAverage, diffSD = diffSD, diffQuantile_25_75 = diffQuantile_25_75, n_indiv = n_indiv,
		n_plots_sampling = ifelse(mode == "spatial", n_plots_sampling, NA)))
}

## Function to recompute indices when subsetting
source("./indices_subsample.R")

#### Load data
## Paths
mainFolder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

clim_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((species_id < 1) | (species_id > length(ls_species)))
	stop(paste0("Species id = ", species_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

species = ls_species[species_id]
print(paste("Script running for species:", species))

treeData = treeData[speciesName_sci == species]
savingPath = paste0("./", species, "/")

if (!dir.exists(savingPath))
	dir.create(savingPath)

## Subsample tree data


# Get parents and children indices
n_indiv = unique(treeData[, .(tree_id, plot_id)])[, .N]
print(paste("Number of individuals:", n_indiv))

if (n_indiv > max_indiv)
{
	print("Too many individuals, subsampling")

	checkSampling = subsampling(treeData, n_indiv_target = max_indiv, mode = "spatial")
	treeData = checkSampling[["sampledData"]]
	n_indiv = checkSampling[["n_indiv"]]

	if (checkSampling[["diffAverage"]])
		warning("The subsample does not look representative of the whole data set, check the average")

	if (checkSampling[["diffSD"]])
		warning("The subsample does not look representative of the whole data set, check the std. dev")

	if (checkSampling[["diffQuantile_25_75"]])
		warning("The subsample does not look representative of the whole data set, check the quantiles 0.25, 0.5, and 0.75")
}

## Climate
climate = readRDS(paste0(clim_folder, "europe_reshaped_climate.rds"))

## Set-up indices
indices = indices_subsample(species_id, run_id, treeData, savingPath, mainFolder, clim_folder)

# indices[plot_id == "504259_france"]
# 504259_france, 1
indices[plot_id == "880396_france"]
# 880396_france, 30

# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]

# Compute the number of growing years per individual
indices[, nbYearsGrowth := max(year) - min(year), by = .(plot_id, tree_id)]

checkUp = all(indices[, nbYearsGrowth == index_clim_end - index_clim_start])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

if (indices[, .N] != treeData[, .N])
	stop(paste0("Dimension mismatch between indices and treeData for species `", species, "`"))

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

nbYearsGrowth = unique(indices[, .(tree_id, plot_id, nbYearsGrowth)])[, nbYearsGrowth]
if (length(nbYearsGrowth) != n_indiv)
	stop("Dimension mismatch between nbYearsGrowth and n_indiv")

n_hiddenState = indices[.N, index_gen]
print(paste("Number of latent dbh:", n_hiddenState))
n_latentGrowth = n_hiddenState - n_indiv
print(paste("Number of latent growth:", n_latentGrowth))

parents_index = treeData[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(plot_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

# if (length(not_parent_index) != indices[.N, index_gen] - n_indiv)
# 	stop("Dimension mismatch between not_parent_index, n_hiddenState, and n_indiv")

#### Compute and save normalising constantes
normalisation(dt = treeData, colnames = "dbh", folder = savingPath, filename = paste0(run_id, "_dbh_normalisation.rds"))
normalisation(dt = climate, colnames = c("pr", "tas"), folder = savingPath, filename = paste0(run_id, "_climate_normalisation.rds"),
	indices = indices, col_ind_start = "index_clim_start", col_ind_end = "index_clim_end")

climate_mu_sd = readRDS(paste0(savingPath, run_id, "_climate_normalisation.rds"))

#### Stan model
## Define stan variables
# Common variables
maxIter = 1500
n_chains = 3

# Initial values for states only
average_yearlyGrowth = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_yearlyGrowth) != n_indiv)
	stop("Dimensions mismatch between average_yearlyGrowth and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
	n_latentGrowth = n_latentGrowth, average_yearlyGrowth = average_yearlyGrowth, nbYearsGrowth = nbYearsGrowth,
	normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]][["latent_growth"]]))

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_climate = climate[, .N], # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_latentGrowth = n_hiddenState - n_indiv, # Dimension of the state space vector for latent growth
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsGrowth = nbYearsGrowth, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	children_index = children_index, # Index of children in the observed data
	climate_index = indices[type == "parent", index_clim_start], # Index of the climate associated to each parent

	# Observations
	Yobs = treeData[, dbh],
	sd_dbh = sd_dbh_beforeSubsample,

	# Explanatory variables
	precip = climate[, pr], # Annual precipitations (sum over 12 months)
	pr_mu = climate_mu_sd[variable == "pr", mu],
	pr_sd = climate_mu_sd[variable == "pr", sd],

	tas = climate[, tas], # Annual average temperature (average over 12 months)
	tas_mu = climate_mu_sd[variable == "tas", mu],
	tas_sd = climate_mu_sd[variable == "tas", sd],

	totalTreeWeight = treeData[parents_index, totalTreeWeight] # Includes also the weight of other species!
)

saveRDS(object = stanData, file = paste0(savingPath, run_id, "_stanData.rds"))

## Compile model
model = cmdstan_model("./growth.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = round(1*maxIter/3), iter_sampling = round(2*maxIter/3), save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 13, adapt_delta = 0.9)

time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-run=", run_id, "-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-run=", run_id, "-", time_ended, ".rds"))

results$cmdstan_diagnose()

results$print(c("lp__", "potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"competition_slope", "measureError", "processError"))


