
#### Aim of prog: plot the log-lik of toyPara_GPUs.stan
## Details:
# I fix all the parameters but the variances to their true values (including the states) and observe the log-lik for 

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
# Set seed
set.seed(1969-08-18) # Woodstock seed

# Number of plots
n_plots = 4

# Number of census per plot
record_min = 3
record_max = 6

n_census_per_plot = sample(x = record_min:record_max, size = n_plots, replace = TRUE)
# n_census_per_plot = rep(3, n_plots)

if (any(n_census_per_plot < 2))
	stop("There should be at least 2 census per plot")

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

sigma_process = 10
sigma_measure = 4

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
	treeStates_dt[index, true_dbh_noError := initial_dbh]

	for (j in 1:(n_census_per_plot[xy] - 1))
	{
		next_dbh = intercept + slope_dbh*treeStates_dt[index + j - 1, true_dbh] +
			slope_precip*treeStates_dt[index + j - 1, precipitations] +
			quad_slope_precip*(treeStates_dt[index + j - 1, precipitations])^2
		treeStates_dt[index + j, true_dbh := rnorm(n = n_indiv_per_plot[xy], mean = next_dbh, sd = sigma_process)]
		treeStates_dt[index + j, true_dbh_noError := next_dbh]
	}
	
	count = count + nb_measures_plot
}

treeStates_dt[, range(true_dbh)]

treeStates_dt[, dbh := rnorm(n = .N, mean = true_dbh, sd = sigma_measure)]
setorder(treeStates_dt, plot_id, tree_id, year)
treeStates_dt[, unique_id := 1:.N]

#### Generate partial data set (keep only first and last measurements)
kept_rows = sort(c(treeStates_dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1],
	treeStates_dt[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]))

treeData = treeStates_dt[kept_rows]
print(paste0(round(treeData[, .N]*100/treeStates_dt[, .N], 2), " % of the data kept"))

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

## Create different variables for indices
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
print(paste("Number of individuals:", n_indiv))

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

Y_generated_0 = rnorm(n_indiv, treeData[parents_index, dbh], 5)



#######################################################?
######## 		THIRD PART: Plot log-lik 		########
#######################################################?

#### Function log lik
loglik_1D = function(variance, list_args)
{
	## Description arguments:
	# 1/ variance is the variable (i.e., 'x'). It is either measureError or processError
	# 2/ list_args is a list containing all the other parameters (set to their true values or not):
	#	- intercepts
	#	- slopes_dbh
	#	- slopes_precip
	#	- quad_slopes_precip
	#	- Y_generated, which are the estimated states by stan
	#	- true_dbh, which are the true dbh from the complete data (i.e., from treeStates_dt[, true_dbh])
	#	- Yobs, which are the observed dbh (i.e., from treeData, incomplete dataset)
	#	- Y_generated_0, which are the initial states for diffuse initialisation
	#	- parents_index, not_parent_index, children_index, which are for the complete dataset
	#	- parentsObs_index, childrenObs_index, which are for the incomplete dataset
	#

	## common variables and check-up
	required_names = c("intercepts", "slopes_dbh", "slopes_precip", "quad_slopes_precip",
		"Y_generated", "true_dbh", "Yobs", "Y_generated_0",
		"parents_index", "not_parent_index", "children_index", "childrenObs_index", "parentsObs_index")
	ls_names = names(list_args)
	matching_names = required_names %in% ls_names

	if (any(!matching_names))
		stop(paste0("The following arguments must be included in list_args:\n- ", paste(required_names[!matching_names], collapse = "\t\n- ")))

	if (!any(c("measureError", "processError") %in% ls_names))
		stop("One variance is required")

	if (all(c("measureError", "processError") %in% ls_names))
	{
		warning("Two variances were provided, but only measureError is kept")
		list_args[["processError"]] = NULL
		ls_names = names(list_args)
	}

	is_measureError_provided = FALSE
	is_processError_provided = FALSE

	if ("measureError" %in% ls_names)
	{
		# print("Plot done for processError")
		is_measureError_provided = TRUE
	}

	if ("processError" %in% ls_names)
	{
		# print("Plot done for measureError")
		is_processError_provided = TRUE
	}

	# The for loop is to vectorised the function, so that it is compatible with 'curve'
	n = length(variance)
	minus_loglik = numeric(n)

	for (i in 1:n)
	{
		target = 0

		## Compute likelihood
		# Priors
		target = target + sum(dnorm(list_args[["intercepts"]], mean = 2, sd = 10, log = TRUE))
		target = target + sum(dgamma(list_args[["slopes_dbh"]], shape = 1.0^2/100, rate = 1.0/100, log = TRUE))
		target = target + sum(dnorm(list_args[["slopes_precip"]], mean = 0, sd = 10, log = TRUE))
		target = target + sum(dnorm(list_args[["quad_slopes_precip"]], mean = 0, sd = 10, log = TRUE))

		if (is_measureError_provided)
		{
			target = target + sum(dgamma(variance[i], shape = 1.0/100, rate = 1.0/100, log = TRUE))
			target = target + sum(dgamma(list_args[["measureError"]], shape = 10.0^2/100, rate = 10.0/100, log = TRUE))
		}

		if (is_processError_provided)
		{
			target = target + sum(dgamma(variance[i], shape = 10.0^2/100, rate = 10.0/100, log = TRUE))
			target = target + sum(dgamma(list_args[["processError"]], shape = 1.0/100, rate = 1.0/100, log = TRUE))
		}

		# Prior on initial hidden state
		target = target + sum(dnorm(list_args[["Y_generated"]][list_args[["parentsObs_index"]]],
			mean = Y_generated_0,
			sd = ifelse(is_processError_provided, list_args[["processError"]], variance[i])), log = TRUE)

		# Observation model
		# Compare true (hidden/latent) parents with observed parents
		target = target + sum(dnorm(list_args[["Yobs"]][parents_index],
			mean = list_args[["Y_generated"]][list_args[["parentsObs_index"]]],
			sd = ifelse(is_measureError_provided, list_args[["measureError"]], variance[i]), log = TRUE))

		# Compare true (hidden/latent) children with observed children
		target = target + sum(dnorm(list_args[["Yobs"]][children_index],
			mean = list_args[["Y_generated"]][list_args[["childrenObs_index"]]],
			sd = ifelse(is_measureError_provided, list_args[["measureError"]], variance[i]), log = TRUE))

		minus_loglik[i] = target;
	}
	return (minus_loglik)
}

#### Plot
## Load the stan results to get the generated states
results = readRDS("./toyPara_GPUs_nonInformative_noVariance_20-08-2021_??h??.rds")
Y_generated_array = results$draws("Y_generated")
class(Y_generated_array)
dim(Y_generated_array)

Y_generated_vector = numeric(dim(Y_generated_array)[3])

for (variable in 1:dim(Y_generated_array)[3])
	Y_generated_vector[variable] = mean(Y_generated_array[,,variable])

## Prepare the arguments for the function loglik_1D
arguments_ls = list(
	# Parameters
	intercepts = intercept,
	slopes_dbh = slope_dbh,
	slopes_precip = slope_precip,
	quad_slopes_precip = quad_slope_precip,

	# Data
	true_dbh = treeStates_dt[, true_dbh],
	Yobs = treeData[, observed_dbh],

	# States
	Y_generated = Y_generated_vector, # treeStates_dt[, true_dbh]
	Y_generated_0 = Y_generated_0,

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	parentsObs_index = indices[type == "parent", index_gen], # Corresponding index of observed parents in Y_generated
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in Y_generated
	not_parent_index = not_parent_index, # Index in Y_generated of states that cannot be compared to data

	# Variances #! Only one variance to provide
	processError = sigma_process
	# measureError = sigma_measure,
	)

loglik_1D(variance = results$draws("measureError")[3000,1,] , list_args = arguments_ls)

curve(loglik_1D(x, list_args = arguments_ls), from = 1.9, to = 2.25,
	xlab = "sigma measure", ylab = "Log-likelihood", lwd = 2)

####? COMMENT AND REFLEXION ZONE

# I verified empirically that the maximum likelihood is not dependent of the diffuse initialisation Y_generated_0

####? END COMMENT AND REFLEXION ZONE

# -------------------------------------------------------

####! CRASH TEST ZONE
plot(treeStates_dt[, true_dbh], Y_generated_vector, xlab = "True dbh",
	ylab = "Estimated dbh", pch = 19, col = "#34568B")

plot(treeData[, dbh], Y_generated_vector[indices[, index_gen]], xlab = "True dbh",
	ylab = "Estimated dbh", pch = 19, col = "#34568B")
abline(a = 0, b = 1, lwd = 2, col = "#CD212A")


test = lm(treeData[, dbh] ~ Y_generated_vector[indices[, index_gen]])

test = lm(Y_generated_vector ~ -1 + treeStates_dt[, true_dbh])

abline(test, lwd = 2, col = "#CD212A")

par(mfrow = c(2, 2))
plot(test)
dev.off()

####! END CRASH TEST ZONE