
# This prog generates data for toyPara_GPUs_1indiv.stan
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var

####
rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(cmdstanr)
library(shiny)

options(max.print = 500)



###########################################################?
######## 		FIRST PART: Create the data 		########
###########################################################?

#### Variables
## Comon variables
n_plots = 1
n_census_per_plot = 418
n_indiv_per_plot = 1

# Number of data
nb_measures = sum(n_indiv_per_plot * n_census_per_plot)

## Storing data
treeStates_dt = data.table(year = integer(nb_measures), tree_id = integer(nb_measures),
	plot_id = integer(nb_measures), true_dbh = numeric(nb_measures))

## Parameters
sigma_process = 10
sigma_measure = 4

#### Generate complete data
init_year_plot = 2005

census_years = init_year_plot:(init_year_plot + n_census_per_plot - 1)
treeStates_dt[, c("year", "plot_id", "tree_id") :=
	.(census_years, 1, 1:n_indiv_per_plot)]

# Initial state
initial_dbh = rgamma(n_indiv_per_plot, shape = 150^2/500, rate = 150/500)
treeStates_dt[1, true_dbh := rnorm(n = 1, mean = initial_dbh, sd = sigma_process)]

for (j in 1:(n_census_per_plot - 1))
{
	next_dbh = treeStates_dt[j, true_dbh]
	treeStates_dt[j + 1, true_dbh := rnorm(n = n_indiv_per_plot, mean = next_dbh, sd = sigma_process)]
}

treeStates_dt[, range(true_dbh)]

treeStates_dt[, observed_dbh := rnorm(n = .N, mean = true_dbh, sd = sigma_measure)]
setorder(treeStates_dt, plot_id, tree_id, year)
treeStates_dt[, unique_id := 1:.N]

#### Generate partial data set (keep random measurements, including first and last)
nb_measures_kept = 50
if (nb_measures_kept <= 2)
{
	warning("Only first and last measurement kept")
	kept_rows = c(1, nb_measures)
}

if (nb_measures_kept > 2)
	kept_rows = sort(c(1, sample(x = 2:(nb_measures - 1), size = nb_measures_kept - 2), nb_measures))

##! Quick and dirty solution to keep everything
# kept_rows = 1:nb_measures

treeData = treeStates_dt[kept_rows]
print(paste0(round(treeData[, .N]*100/treeStates_dt[, .N], 2), " % of the data kept"))

if (isTRUE(all.equal(treeData, treeStates_dt)))
	warning("There is no missing data, are you sure it is what you want?")

setorder(treeData, plot_id, tree_id, year)

#### plots
## General plot
plot(treeStates_dt[, unique_id], treeStates_dt[, true_dbh], pch = 19, cex = 0.7, col = "red", ty = "o", xlab = "unique id", ylab = expression(observed_dbh[i](t)),
	ylim = c(min(treeStates_dt[, .(true_dbh, observed_dbh)]) - abs(min(treeStates_dt[, .(true_dbh, observed_dbh)]))/10,
		max(treeStates_dt[, .(true_dbh, observed_dbh)] + max(treeStates_dt[, .(true_dbh, observed_dbh)]/10))), las = 1)
points(treeData[, unique_id], treeData[, observed_dbh], pch = 3, cex = 0.8, col = "blue", ty = "o", lty = 3)

legend("top", legend = c("Obs. (with missing data)", "True states"), pch = c(3, 19),
	col = c("blue", "red"), lty = c(3, 1), horiz=TRUE, bty="n", cex=0.9)



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

## Data table indices
indices = data.table(year = integer(nb_measures), tree_id = integer(nb_measures),
	plot_id = integer(nb_measures), index_gen = integer(nb_measures))

years_indices = fillYears(treeData[, year])
indices[, year := years_indices[["fill_years"]]]
indices[, tree_id := 1]
indices[, plot_id := 1]
indices[years_indices[["indices"]], index_gen := years_indices[["indices"]]]

indices = indices[index_gen != 0]

## Set everyone to child type
indices[, type := "child"]

## Correct for those who are parent type (Should be the first line only!)
indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]

## Compute nb years per individual
indices[, nbYearsPerIndiv := max(year) - min(year) + 1, by = .(plot_id, tree_id)]

if (indices[, .N] != treeData[, .N])
	stop("Dimensions mismatch between indices and treeData")

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

n_indiv = unique(treeData[, .(tree_id, plot_id)])[, .N]
print(paste("Number of individuals:", n_indiv))

nbYearsPerIndiv = unique(indices[, .(tree_id, plot_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
if (length(nbYearsPerIndiv) != n_indiv)
	stop("Dimension mismatch between nbYearsPerIndiv and n_indiv")

## Create the required indices vectors used in stan
parents_index = treeData[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(plot_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]
not_parent_index = 1:indices[.N, index_gen]
not_parent_index = not_parent_index[!(not_parent_index %in% indices[type == "parent", index_gen])]

if (length(parents_index) != n_indiv)
	stop("Dimensions mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimensions mismatch between children_index and number of children")

if (length(not_parent_index) != indices[.N, index_gen] - n_indiv)
	stop("Dimensions mismatch between not_parent_index, n_hiddenState, and n_indiv")



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
		Y_gen[count + 1] = rnorm(1, mean = dbh_parents[i], sd = 5)
		for (j in 2:years_indiv[i])
			Y_gen[count + j] = Y_gen[count + j - 1] + rnorm(1, mean = 0, sd = 10)

		count = count + years_indiv[i];
	}

	return(list(Y_generated = Y_gen))
}

Y_generated_0 = rnorm(n_indiv, treeData[parents_index, observed_dbh], 5)

#### Stan model
## Define stan variables
# Common variables
maxIter = 4e3
n_chains = 3

range(parents_index)
range(indices[type == "parent", index_gen])  
range(children_index)
range(indices[type == "child", index_gen])
range(not_parent_index)

stanData = list(
	# Number of data
	n_obs = n_obs, # Number of trees observations
	n_hiddenState = indices[.N, index_gen], # Dimension of the state space vector
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsPerIndiv = nbYearsPerIndiv, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	parentsObs_index = indices[type == "parent", index_gen], # Corresponding index of observed parents in Y_generated
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in Y_generated
	not_parent_index = not_parent_index, # Index in Y_generated of states that cannot be compared to data

	# Initial state, used to start the states and compare with Y_generated[parentsObs_index] (i.e., a prior parameter)
	Y_generated_0 = Y_generated_0,

	# Observations
	Yobs = treeData[, observed_dbh]
)

# Initial value for states only
average_G = (treeData[last_child_index, observed_dbh] - treeData[parents_index, observed_dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, observed_dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen])

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

## Compile stan model
model = cmdstan_model("toyPara_GPUs_1indiv.stan")

## Run model
# start = proc.time()

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains, # threads_per_chain = 2,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, init = initVal_Y_gen, save_warmup = TRUE,
	max_treedepth = 14, adapt_delta = 0.95)

# proc.time() - start

# results$save_object(file = paste0("./toyPara_GPUs_nonInformative_noProcessError_", format(Sys.time(), "%d-%m-%Y_%Hh%M"), ".rds"))

results$cmdstan_diagnose()

#### Plots of posterior distributions and traces
# ## Intercept
# plot_title = ggplot2::ggtitle("Traces for intercept")
# mcmc_trace(results$draws("intercepts")) + plot_title

# plot_title = ggplot2::ggtitle("Posterior distribution intercept", "with medians and 80% intervals")
# mcmc_areas(results$draws("intercepts"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = intercept, color = "#FFA500")

# ## Slope for dbh
# plot_title = ggplot2::ggtitle("Posterior distribution slope dbh", "with medians and 80% intervals")
# mcmc_areas(results$draws("slopes_dbh"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = slope_dbh, color = "#FFA500")

# plot_title = ggplot2::ggtitle("Traces for slope dbh")
# mcmc_trace(results$draws("slopes_dbh")) + plot_title

# ## Slope for precipitation
# plot_title = ggplot2::ggtitle("Posterior distribution slope precip", "with medians and 80% intervals")
# mcmc_areas(results$draws("slopes_precip"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = slope_precip, color = "#FFA500")

# plot_title = ggplot2::ggtitle("Traces for slope precip")
# mcmc_trace(results$draws("slopes_precip")) + plot_title

# ## Slope for precipitation (quadratique term)
# plot_title = ggplot2::ggtitle("Posterior distribution slope precip (quadratic term)", "with medians and 80% intervals")
# mcmc_areas(results$draws("quad_slopes_precip"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = quad_slope_precip, color = "#FFA500")

# plot_title = ggplot2::ggtitle("Traces for slope precip (quadratic term)")
# mcmc_trace(results$draws("quad_slopes_precip")) + plot_title

## Measurement error
plot_title = ggplot2::ggtitle("Posterior distribution measure error", "with medians and 80% intervals")
mcmc_areas(results$draws("measureError"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigma_measure, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for measure error")
mcmc_trace(results$draws("measureError")) + plot_title

## Process error
plot_title = ggplot2::ggtitle("Posterior distribution process error", "with medians and 80% intervals")
mcmc_areas(results$draws("processError"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigma_process, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for process error")
mcmc_trace(results$draws("processError", inc_warmup = TRUE)) + plot_title

# intercept = 2.4

# slope_precip = 0.004
# quad_slope_precip = -0.00001

# slope_dbh = 1.1

# sigma_process = 4
# sigma_measure = 1

## Check-up states
chosen_state = 50
if (chosen_state > treeStates_dt[, .N])
	stop(paste("chosen_state must be smaller than", treeStates_dt[, .N]))

plot_title = ggplot2::ggtitle(paste("Traces for state", chosen_state))
mcmc_trace(results$draws(paste0("Y_generated[", chosen_state, "]"))) + plot_title

plot_title = ggplot2::ggtitle(paste("Posterior distribution for state", chosen_state), "with medians and 80% intervals")
mcmc_areas(results$draws(paste0("Y_generated[", chosen_state, "]")), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = treeStates_dt[chosen_state, true_dbh], color = "#FFA500")
