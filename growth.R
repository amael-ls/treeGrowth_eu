
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
# library(doParallel)
library(cmdstanr)
library(callr)

# #### Create the cluster
# ## Cluster variables
# nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

# print("nodeslist")
# print(nodeslist)
# print("end nodeslist")

# ## Make cluster
# cl = makeCluster(nodeslist, type = "PSOCK")
# registerDoParallel(cl)

# print("cluster done")

# print(paste("number of cores:", future::availableCores()))

#### Tool function
## Initiate Y_gen with reasonable value (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "years_indiv", "average_G", "n_hiddenState")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide chain_id, Yobs, and years_indiv")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	years_indiv = providedArgs[["years_indiv"]]
	average_G = providedArgs[["average_G"]]
	n_hiddenState = providedArgs[["n_hiddenState"]]

	Y_gen = numeric(n_hiddenState)

	count = 0

	for (i in 1:n_indiv) # Not that this forbid trees to shrink
	{
		Y_gen[count + 1] = rgamma(1, shape = dbh_parents[i]^2, rate = dbh_parents[i]) # Mean = dbh_parents[i], Variance = 1
		for (j in 2:years_indiv[i])
			Y_gen[count + j] = Y_gen[count + j - 1] + average_G[i] + rgamma(1, shape = 0.25, rate = 0.5)

		count = count + years_indiv[i];
	}

	return(list(Y_generated = Y_gen))
}

#### Data
if (!file.exists("./tilia_growthData.rds"))
	script(script = "./create_tilia_growthData.R")

treeData = readRDS("./tilia_growthData.rds")


# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]

# Compute nb years per individual
indices[, nbYearsPerIndiv := max(year) - min(year) + 1, by = .(plot_id, tree_id)]

checkUp = all(indices[, nbYearsPerIndiv == index_precip_end - index_precip_start + 1])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

if (length(precip) != indices[.N, index_precip_end])
	stop("Dimension mismatch between climate and indices")

if (indices[, .N] != treeData[, .N])
	stop("Dimension mismatch between indices and treeData")

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

n_indiv = unique(treeData[, .(tree_id, pointInventory_id)])[, .N]
print(paste("Number of individuals:", treeData[, .N]))

nbYearsPerIndiv = unique(indices[, .(tree_id, plot_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
if (length(nbYearsPerIndiv) != n_indiv)
	stop("Dimension mismatch between nbYearsPerIndiv and n_indiv")

parents_index = treeData[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(pointInventory_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(pointInventory_id, tree_id)][, V1]
not_parent_index = 1:indices[.N, index_gen]
not_parent_index = not_parent_index[!(not_parent_index %in% indices[type == "parent", index_gen])]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

if (length(not_parent_index) != indices[.N, index_gen] - n_indiv)
	stop("Dimension mismatch between not_parent_index, n_hiddenState, and n_indiv")

#### Stan model
## Define stan variables
# Common variables
maxIter = 2e2 # 2e3
n_chains = 1 # 4

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_precip = length(precip), # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_hiddenState = indices[.N, index_gen], # Dimension of the state space vector
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsPerIndiv = nbYearsPerIndiv, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	parentsObs_index = indices[type == "parent", index_gen], # Corresponding index of observed parents in Y_generated
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in Y_generated
	climate_index = indices[type == "parent", index_precip_start], # Index of the climate associated to each parent
	not_parent_index = not_parent_index, # Index in Y_generated of states that cannot be compared to data

	# Observations
	Yobs = treeData[, dbh],

	# Explanatory variable
	precip = precip, # Precipitations

	# Parameter for parralel calculus
	grainsize = 1
)

stanfit <- rstan::read_stan_csv(fit$output_files())