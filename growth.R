
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)
#! IMPORTANT REMARK:
# There are two ways to use Y_generated_0:
#	1. It is a data, and therefore can be used at any time to help convergence on either the parents only or on the hidden states too
#		depending on the size of the vector and the likelihood
#	2. It is an initial condition, which help to start the hidden states. This also means that Y_generated_0 "is forgotten" as iterations
#		goes by
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(doParallel) #! TO UNCOMMENT FOR COMPUTE CANADA
library(bayesplot)
library(cmdstanr)
library(stringi)

if (!("callr" %in% installed.packages()))
	install.packages("callr", repos = "http://cran.us.r-project.org")

if (!("future" %in% installed.packages()))
	install.packages("future", repos = "http://cran.us.r-project.org")

#### Create the cluster
## Cluster variables
# array_id = 150 #! REMOVE THIS LINE WHEN DONE WITH TEST
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) # Each array takes care of one year and process all the variables for that year
print(paste0("array id = ", array_id))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

print(paste("number of cores:", future::availableCores()))

#### Tool function
## Initiate Y_gen with reasonable value (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "years_indiv", "average_G", "n_hiddenState", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide chain_id, Yobs, years_indiv, and normalise")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	years_indiv = providedArgs[["years_indiv"]]
	average_G = providedArgs[["average_G"]]
	n_hiddenState = providedArgs[["n_hiddenState"]]
	normalise = providedArgs[["normalise"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	Y_gen = numeric(n_hiddenState)

	count = 0

	for (i in 1:n_indiv) # Not that this forbid trees to shrink
	{
		Y_gen[count + 1] = rgamma(1, shape = dbh_parents[i]^2, rate = dbh_parents[i]) # Mean = dbh_parents[i], Variance = 1
		for (j in 2:years_indiv[i])
			Y_gen[count + j] = Y_gen[count + j - 1] + average_G[i] + rgamma(1, shape = 2.5, rate = 5) # mean = 0.5, var = 0.1

		count = count + years_indiv[i];
	}

	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
	}
	return(list(latentState = Y_gen))
}

## Function to ... #! works exclusively when trees are measured twice (like in the French data)
forced_states = function(dbh_parents, dbh_children, parentsObs_index, childrenObs_index, years_indiv, n_hiddenState)
{
	states = numeric(n_hiddenState)
	other_index = 1:n_hiddenState
	other_index = other_index[!(other_index %in% parentsObs_index) & !(other_index %in% childrenObs_index)]
	
	interp_slop = (dbh_children - dbh_parents)/years_indiv
	
	states[parentsObs_index] = dbh_parents
	states[childrenObs_index] = dbh_children

	count = 0
	for (indiv in 1:length(dbh_parents))
	{
		for (yr in 1:(years_indiv[indiv] - 1))
		{
			count = count + 1
			states[other_index[count]] = dbh_parents + yr*interp_slop[indiv]
		}
	}

	return(states)
}

#### Load data
## Paths
mainFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

clim_folder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/climateData/Chelsa/yearlyAverage/"

if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((array_id < 1) | (array_id > length(ls_species)))
	stop(paste0("Array id = ", array_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

species = ls_species[array_id]
treeData = treeData[speciesName_sci == species]
savingPath = paste0("./", species, "/")
if (!dir.exists(savingPath))
	dir.create(savingPath)

## Climate
climate = readRDS(paste0(clim_folder, "FR_reshaped_climate.rds"))

## inidices
indices = readRDS(paste0(mainFolder, species, "_indices.rds"))

# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1], type := "parent"]

# Compute nb years per individual
indices[, nbYearsPerIndiv := max(year) - min(year) + 1, by = .(pointInventory_id, tree_id)]

checkUp = all(indices[, nbYearsPerIndiv == index_clim_end - index_clim_start + 1])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

if (climate[, .N] != indices[.N, index_clim_end]) #! THIS LINE IS NOW WRONG I THINK! INDEED, CLIMATE IS NOT SPECIES SPECIFIC!
	stop("Dimension mismatch between climate and indices")

if (indices[, .N] != treeData[, .N])
	stop(paste0("Dimension mismatch between indices and treeData for species `", species, "`"))

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

n_indiv = unique(treeData[, .(tree_id, pointInventory_id)])[, .N]
print(paste("Number of individuals:", n_indiv))

nbYearsPerIndiv = unique(indices[, .(tree_id, pointInventory_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
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
maxIter = 3e3
n_chains = 3

# Initial value for states only
average_G = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen], normalise = TRUE, mu_dbh = mean(treeData[, dbh]), sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

# jpeg("./initVal.jpg", height = 1080, width = 1080, quality = 100)
# plot(1:indices[.N, index_gen], initVal_Y_gen[[1]]$latentState, pch = 19, col = "#34568B",
# 	xlab = "Tree index", ylab = "Diameter at breast height (in mm)")
# points(x = indices[type == "parent", index_gen], y = treeData[parents_index, dbh], pch = 19, col = "#FA7A35")
# points(x = indices[type == "child", index_gen], y = treeData[children_index, dbh], pch = 19, col = "#CD212A")
# dev.off()

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_precip = climate[, .N], # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_hiddenState = indices[.N, index_gen], # Dimension of the state space vector
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsPerIndiv = nbYearsPerIndiv, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	parentsObs_index = indices[type == "parent", index_gen], # Corresponding index of observed parents in latentState
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in latentState
	climate_index = indices[type == "parent", index_clim_start], # Index of the climate associated to each parent
	not_parent_index = not_parent_index, # Index in latentState of states that cannot be compared to data

	# Observations
	Yobs = treeData[, dbh],

	# Explanatory variable
	precip = climate[, pr], # Precipitations

	# Diffuse initialisation for the parents #! Need to think more about it, not used for now
	# Y_generated_0 = rnorm(n = n_indiv, treeData[parents_index, dbh], sd = 1),
	initialParents = initVal_Y_gen[[1]]$latentState[indices[type == "parent", index_gen]]

	# Parameter for parallel calculus
	# grainsize = 1
)

## Compile model
model = cmdstan_model("growth.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 150, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 12, adapt_delta = 0.8)

time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-", format(Sys.time(), "%Y-%m-%d_%Hh%M"), ".rds"))

diagnose = results$cmdstan_diagnose()

# stanfit <- rstan::read_stan_csv(fit$output_files())
#### Few plots
# Folder to save plots
figurePath = paste0(savingPath, time_ended, "/")
if (!dir.exists(figurePath))
	dir.create(figurePath)

# Intercept
plot_title = ggplot2::ggtitle("Posterior distribution intercepts", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("intercepts"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "intercept-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for intercepts")
temp_plot = mcmc_trace(results$draws("intercepts")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "intercept-traces.pdf"), plot = temp_plot, device = "pdf")

# Slopes_dbh
plot_title = ggplot2::ggtitle("Posterior distribution slopes_dbh", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("slopes_dbh"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "slopes_dbh-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for slopes_dbh")
temp_plot = mcmc_trace(results$draws("slopes_dbh")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "slopes_dbh-traces.pdf"), plot = temp_plot, device = "pdf")

# Slopes_precip
plot_title = ggplot2::ggtitle("Posterior distribution slopes_precip", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("slopes_precip"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "slopes_precip-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for slopes_precip")
temp_plot = mcmc_trace(results$draws("slopes_precip")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "slopes_precip-traces.pdf"), plot = temp_plot, device = "pdf")

# Quad_slopes_precip
plot_title = ggplot2::ggtitle("Posterior distribution quad_slopes_precip", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("quad_slopes_precip"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "quad_slopes_precip-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for quad_slopes_precip")
temp_plot = mcmc_trace(results$draws("quad_slopes_precip")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "quad_slopes_precip-traces.pdf"), plot = temp_plot, device = "pdf")

# processError
plot_title = ggplot2::ggtitle("Posterior distribution processError", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("processError"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "processError-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for processError")
temp_plot = mcmc_trace(results$draws("processError")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "processError-traces.pdf"), plot = temp_plot, device = "pdf")

# measureError
plot_title = ggplot2::ggtitle("Posterior distribution measureError", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("measureError"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "measureError-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for measureError")
temp_plot = mcmc_trace(results$draws("measureError")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "measureError-traces.pdf"), plot = temp_plot, device = "pdf")

for (i in c(1:3, 12:14))
{
	plot_title = ggplot2::ggtitle(paste0("Posterior distribution latentState[", i, "]"),
		"with medians and 80% intervals")
	temp_plot = mcmc_areas(results$draws(paste0("latentState[", i, "]")), prob = 0.8) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latentState_", i, "-distrib.pdf"), plot = temp_plot, device = "pdf")

	plot_title = ggplot2::ggtitle(paste0("Traces for latentState[", i, "]"))
	temp_plot = mcmc_trace(results$draws(paste0("latentState[", i, "]"), inc_warmup = TRUE)) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latentState_", i, "-traces.pdf"), plot = temp_plot, device = "pdf")
}

####! CRASH TEST ZONE
## Create data
# n = 1e4
# intercept = 2
# slope = 0.5
# sigma = 2.3

# x = runif(n, 0, 100)
# x_norm = (x - mean(x))/sd(x)

# y = rnorm(n, intercept + slope*x, sigma)
# y_norm = (y - mean(y))/sd(y)

# ## Run linear model, no normalisation
# lin = lm(y ~ x)
# mean(residuals(lin)) # Should be near 0
# sd(residuals(lin))/sigma # Should be near 1
# coefficients(lin)["(Intercept)"]/intercept # Should be near 1
# coefficients(lin)["x"]/slope # Should be near 1

# ## Run linear model, x is normalised
# lin_norm = lm(y ~ x_norm)
# mean(residuals(lin_norm)) # Should be near 0
# sd(residuals(lin_norm))/sigma # Should be near 1
# beta_0 = coefficients(lin_norm)["(Intercept)"]
# beta_1 = coefficients(lin_norm)["x_norm"]
# (beta_0 - beta_1*mean(x)/sd(x))/intercept # Should be near 1
# beta_1/sd(x)/slope # Should be near 1

# ## Run linear model, y is normalised
# lin_norm2 = lm(y_norm ~ x)
# mean(residuals(lin_norm2)) # Should be near 0
# sd(residuals(lin_norm2))*sd(y)/sigma # Should be near 1
# beta_0 = coefficients(lin_norm2)["(Intercept)"]
# beta_1 = coefficients(lin_norm2)["x"]
# (sd(y)*beta_0 + mean(y))/intercept # Should be near 1
# beta_1*sd(y)/slope # Should be near 1

# ## Run linear model, both x and y are normalised
# lin_norm3 = lm(y_norm ~ x_norm)
# mean(residuals(lin_norm2)) # Should be near 0
# sd(residuals(lin_norm2))*sd(y)/sigma # Should be near 1
# beta_0 = coefficients(lin_norm3)["(Intercept)"]
# beta_1 = coefficients(lin_norm3)["x_norm"]
# (sd(y)*beta_0 + mean(y) - mean(x)/sd(x)*sd(y)*beta_1)/intercept # Should be near 1
# (sd(y)*beta_1)/sd(x)/slope # Should be near 1

####! END CRASH TEST ZONE
