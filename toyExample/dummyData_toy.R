
# This prog generates data for toy.stan, toyPara_GPUs.stan and toyPara_reduce_sum.stan
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var

####
rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(cmdstanr)

options(max.print = 500)
###########################################################?
######## 		FIRST PART: Create the data 		########
###########################################################?

#### Variables
## Comon variables
# Set seed
set.seed(1969-08-18) # Woodstock seed

# Number of plots
n_plots = 4

# Number of census per plot
record_min = 4
record_max = 8

n_census_per_plot = sample(x = record_min:record_max, size = n_plots, replace = TRUE)

# Number of individuals per plot
indiv_min = 10
indiv_max = 30
n_indiv_per_plot = sample(x = indiv_min:indiv_max, size = n_plots, replace = TRUE)

# Number of data
nb_measures = sum(n_indiv_per_plot * n_census_per_plot)

## Storing data
treeStates_dt = data.table(year = integer(nb_measures), tree_id = integer(nb_measures),
	plot_id = integer(nb_measures), true_dbh = numeric(nb_measures), precipitations = numeric(nb_measures))

## Parameters
intercept = 2.4

slope_precip = 0.004
quad_slope_precip = -0.00001

slope_dbh = 1.1

sigma_process = 4
sigma_measure = 1

#### Generate complete data
init_year_plot = sample(1990:2005, size = n_plots, replace = TRUE)

count = 0
for (xy in 1:n_plots)
{
	census_years = rep(init_year_plot[xy]:(init_year_plot[xy] + n_census_per_plot[xy] - 1), n_indiv_per_plot[xy])
	nb_measures_plot = n_census_per_plot[xy] * n_indiv_per_plot[xy]
	treeStates_dt[(count + 1):(count + nb_measures_plot), c("year", "plot_id", "tree_id") :=
		.(census_years, xy, rep(1:n_indiv_per_plot[xy], each = n_census_per_plot[xy]))]

	precip = runif(n = n_census_per_plot[xy], min = 650, max = 1200)

	treeStates_dt[(count + 1):(count + nb_measures_plot), precipitations := rep(precip, n_indiv_per_plot[xy])]

	# Initial state
	initial_dbh = rgamma(n_indiv_per_plot[xy], shape = 150^2/300, rate = 150/300)
	index = seq(from = 1, to = nb_measures_plot, by = n_census_per_plot[xy]) + count
	treeStates_dt[index, true_dbh := initial_dbh]

	for (j in 1:(n_census_per_plot[xy] - 1))
	{
		next_dbh = intercept + slope_dbh*treeStates_dt[index + j - 1, true_dbh] +
			slope_precip*treeStates_dt[index + j - 1, precipitations] +
			quad_slope_precip*(treeStates_dt[index + j - 1, precipitations])^2
		treeStates_dt[index + j, true_dbh := rnorm(n = n_indiv_per_plot[xy], mean = next_dbh, sd = sigma_process)]
	}
	
	count = count + nb_measures_plot
}

treeStates_dt[, range(true_dbh)]

treeStates_dt[, dbh := rnorm(n = .N, mean = true_dbh, sd = sigma_measure)]

#### Generate partial data set (keep only first and last measurements)
kept_rows = sort(c(treeStates_dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1],
	treeStates_dt[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]))
treeData = treeStates_dt[kept_rows]
setorder(treeData, plot_id, tree_id, year)



###########################################################?
######## 		Second PART: Create indices 		########
###########################################################?

#### Tool function
fillYears = function(years)
{
	if (length(years) < 2)
		stop("From fillYears: Their should be at least two years to fill the gaps")

	if (is.unsorted(years))
		stop("From fillYears: years are assumed to be sorted!")

	fill_years = years[1]:years[length(years)]
	indices = which(fill_years %in% years)
	
	return (list(fill_years = fill_years, indices = indices))
}

#### Get precipitations (in real, from a raster but here from treeStates_dt)
values = unique(treeStates_dt[, .(plot_id, year, precipitations)])
setnames(values, old = "plot_id", new = "id")

#### 'Joining' climate with tree data
count = 0
start = 0
end = 0
iter = 0

nbIndiv = unique(treeData[, .(plot_id, tree_id)])[, .N]
length_filled_years = sum(treeData[, max(year) - min(year) + 1, by = .(plot_id, tree_id)][, V1])

indices = data.table(year = integer(length_filled_years), tree_id = integer(length_filled_years),
	plot_id = integer(length_filled_years), index_gen = integer(length_filled_years),
	index_precip_start = integer(length_filled_years), index_precip_end = integer(length_filled_years))

for (plot in treeData[, unique(plot_id)])
{
	for (indiv in treeData[plot_id == plot, unique(tree_id)])
	{
		years_indices = fillYears(treeData[plot_id == plot & tree_id == indiv, year])
		start = end + 1
		end = end + length(years_indices[["fill_years"]])
		indices[start:end, year := years_indices[["fill_years"]]]
		indices[start:end, tree_id := indiv]
		indices[start:end, plot_id := plot]
		indices[years_indices[["indices"]] + count, index_gen := years_indices[["indices"]] + count]
		count = count + years_indices[["indices"]][length(years_indices[["indices"]])]
		iter = iter + 1
		if (iter %% 1000 == 0)
			print(paste(round(iter*100/nbIndiv, digits = 3), "% done"))
	}
}

## Create the indices for the climate, and format climate data
start = 0
end = 0
count = 0
iter = 0

length_clim = sum(treeData[, max(year) - min(year) + 1, by = .(plot_id)][, V1])
precipitations_yearly = numeric(length_clim)

for (plot in indices[, unique(plot_id)])
{
	precip_years = indices[plot_id == plot, sort(unique(year))]
	
	start = end + 1
	end = start + length(precip_years) - 1

	precipitations_yearly[start:end] = as.numeric(values[id == plot & year %in% precip_years , precipitations])

	for (tree in indices[plot_id == plot, unique(tree_id)])
	{
		precip_start = min(which(indices[tree_id == tree & plot_id == plot, year] %in% precip_years)) + count
		precip_end = max(which(indices[tree_id == tree & plot_id == plot, year] %in% precip_years)) + count

		indices[tree_id == tree & plot_id == plot,
			c("index_precip_start", "index_precip_end") := .(precip_start, precip_end)]
		
		iter = iter + 1
		if (iter %% 1000 == 0)
			print(paste(round(iter*100/nbIndiv, digits = 3), "% done"))
	}
	count = count + length(precip_years)
}

indices = indices[index_gen != 0]



###########################################################?
######## 		THIRD PART: Run stan program 		########
###########################################################?

#### Tool function
## Initiate Y_gen with reasonable value (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "years_indiv", "average_G", "n_hiddenState")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide chain_id, Yobs, and years_indiv")
	
	# print(providedArgs)

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

#### Load data
## Data are already in memory
precip = precipitations_yearly
rm(precipitations_yearly)

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

n_indiv = unique(treeData[, .(tree_id, plot_id)])[, .N]
print(paste("Number of individuals:", treeData[, .N]))

nbYearsPerIndiv = unique(indices[, .(tree_id, plot_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
if (length(nbYearsPerIndiv) != n_indiv)
	stop("Dimension mismatch between nbYearsPerIndiv and n_indiv")

parents_index = treeData[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(plot_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]
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
maxIter = 1e4
n_chains = 4

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

# Initial value for states only
average_G = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen])

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

## Compile stan model
# model = stan_model(file = "./toy.stan")
# model = cmdstan_model("toyPara_reduce_sum.stan", cpp_options = list(stan_threads = TRUE)) # list(stan_threads = TRUE, stan_opencl = TRUE)
# model = cmdstan_model("toyPara_GPUs.stan") #, cpp_options = list(stan_opencl = TRUE))

model = cmdstan_model("toyPara_GPUs.stan", cpp_options = list(stan_threads = TRUE))

## Run model
start = proc.time()

# results = stan(file = "toy.stan", data = stanData, cores = n_chains,
# 	iter = maxIter, chains = n_chains, init = initVal_Y_gen)

# results = model$sample(data = stanData, parallel_chains = n_chains, threads_per_chain = 8, refresh = 2,
# 	iter_warmup = maxIter/2, iter_sampling = maxIter/2, chains = n_chains, init = initVal_Y_gen)

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 1000, chains = n_chains,
	threads_per_chain = 2, iter_warmup = maxIter/2, iter_sampling = maxIter/2, init = initVal_Y_gen)

# results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 2, chains = n_chains,
# 	opencl_ids = c(0, 0), iter_warmup = maxIter/2, iter_sampling = maxIter/2, init = initVal_Y_gen)

proc.time() - start

results

plot_title = ggplot2::ggtitle("Posterior distributions", "with medians and 80% intervals")

# "slopes_dbh", "slopes_precip", "quad_slopes_precip", "processError", "measureError"
mcmc_areas(results$draws("intercepts")) + plot_title

# intercept = 2.4

# slope_precip = 0.004
# quad_slope_precip = -0.00001

# slope_dbh = 1.1

# sigma_process = 4
# sigma_measure = 1