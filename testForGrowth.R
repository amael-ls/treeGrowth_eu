
#### Aim of prog: To test the model growth_2ndOption.stan with dummy data

# https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
library(DHARMa)

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
			Y_gen[count + j] = Y_gen[count + j - 1] + rgamma(1, shape = average_G[i]^2/5, rate = average_G[i]/5) # mean = 0.5, var = 0.1

		count = count + years_indiv[i];
	}

	Y_gen = rgamma(n_indiv, dbh_parents^2, dbh_parents)
	
	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
	}
	# return(list(latent_dbh = Y_gen, potentialGrowth = 0, dbh_slope = 0, pr_slope = 0, pr_slope2 = 0))
	return(list(latent_dbh_parents = Y_gen))
}

## Log-likelihood function for test 4
loglik = function(data, latent_dbh_parents, latent_dbh_children, parents_index, children_index, nbYearsPerIndiv, params, normalised_precip,
	climate_index, measureError, processError)
{
	n_indiv = length(parents_index)
	computed_latent_child = numeric(length = n_indiv)

	potentialGrowth = params[["potentialGrowth"]]
	dbh_slope = params[["dbh_slope"]]
	pr_slope = params[["pr_slope"]]
	pr_slope2 = params[["pr_slope2"]]

	for (i in 1:n_indiv)
	{
		next_dbh = latent_dbh_parents[i];
		for (j in 2:nbYearsPerIndiv[i])
			next_dbh = next_dbh + growth_fct(next_dbh, normalised_precip[climate_index[i] + j - 2], params);
		
		computed_latent_child[i] = next_dbh;
	}
	
	ll = sum(dnorm(data[parents_index, dbh], mean = latent_dbh_parents, sd = measureError, log = TRUE)) +
		sum(dnorm(data[children_index, dbh], mean = computed_latent_child, sd = measureError, log = TRUE)) +
		dnorm(params["potentialGrowth"], mean = 0, sd = 100, log = TRUE) + dnorm(params["dbh_slope"], mean = 0, sd = 5, log = TRUE) +
		sum(dnorm(pr_slope, mean = 0, sd = 5, log = TRUE)) + sum(dnorm(pr_slope2, mean = 0, sd = 5, log = TRUE)) +
		sum(dnorm(measureError, mean = 3, sd = 0.25, log = TRUE))
	return (ll)
}

## Function to compute the residuals (latent state)
myPredictObs = function(draws_array, id_latent, regex = "latent_dbh")
{
	id_latent = unique(id_latent)
	if (length(id_latent) != 1)
		stop("A single id should be provided")
	
	n_chains = ncol(draws_array)
	length_chain = nrow(draws_array)
	output = numeric(length = n_chains*length_chain)

	for (i in 1:n_chains)
	{
		start = (i - 1)*length_chain + 1
		end = i*length_chain
		output[start:end] = draws_array[, i, paste0(regex, "[", id_latent, "]")]
	}
	return (output)
}

## Compute next dbh
growth_fct = function(dbh0, precip, params)
{
	potentialGrowth = params["potentialGrowth"]
	dbh_slope = params["dbh_slope"]
	pr_slope = params["pr_slope"]
	pr_slope2 = params["pr_slope2"]
	return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2))
}

## Bayesplot is having troubles on my mac (Arial font not always found), so I create my own traces plot
lazyTrace = function(draws, ...)
{
	if (!all(class(draws) %in% c("draws_array", "draws", "array")))
		stop("The class of draws should be compatible with stan (draws_array, draws, array)")
	
	n_chains = dim(draws)[2]
	n_iter = dim(draws)[1]
	cols = c("#845D29", "#D8C29D", "#178F92", "#1D1F54")

	min_val = min(draws)
	max_val = max(draws)

	providedArgs = list(...)
	nbArgs = length(providedArgs)

	ls_names = names(providedArgs)

	val_ind = stri_detect(str = ls_names, regex = "val")
	xlab_ind = (ls_names == "xlab")
	ylab_ind = (ls_names == "ylab")
	main_ind = (ls_names == "main")

	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	plot(0, pch = "", xlim = c(0, n_iter), ylim = scaling*c(min_val, max_val), axes = TRUE, bg = "transparent",
		xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""),
		ylab = ifelse(any(ylab_ind), providedArgs[["ylab"]], ""),
		main = ifelse(any(main_ind), providedArgs[["main"]], ""))

	for (chain in 1:n_chains)
		lines(1:n_iter, scaling*draws[, chain,], type = "l", col = cols[chain])
		# lines(1:n_iter, draws[, chain,], type = "l", col = cols[chain])
	if (any(val_ind))
		{
			for (val in ls_names[val_ind])
				abline(h = scaling*providedArgs[[val]], col = "#CD212A", lwd = 4)
		}
}

#### Create dummy data
## Load real data
treeData = readRDS("Tilia_platyphyllos/trees_forest_reshaped.rds")
species = "Tilia_platyphyllos"
treeData = treeData[speciesName_sci == species]
climate = readRDS("Tilia_platyphyllos/FR_reshaped_climate.rds")
indices = readRDS("Tilia_platyphyllos/Tilia_platyphyllos_indices.rds")

climate_mu_sd = readRDS("Tilia_platyphyllos/climate_normalisation.rds")

mu_pr = climate_mu_sd[variable == "pr", mu]
sd_pr = climate_mu_sd[variable == "pr", sd]

## Set-up indices
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

## Parameters (to be estimated). /!\ Given on a non standardised scale for dbh
potentialGrowth = log(2.47); # i.e., in average a tree grows 2.47 cm
dbh_slope = 0.11/135;

pr_slope = -0.34;
pr_slope2 = 0.088;

params = c(potentialGrowth = potentialGrowth, dbh_slope = dbh_slope, pr_slope = pr_slope, pr_slope2 = pr_slope2)

processError = 5.68 # This is a variance: var /!\
measurementError = sqrt(1) # This is a standard deviation: sd /!\ i.e., the variance is 3

## Simulate dbh
set.seed(123)
n = treeData[, .N]/2
if (n != length(children_index))
	stop("Length mismatch between n and children_index")

latent_dbh = data.table(dbh0 = treeData[parents_index, dbh], dbh1 = numeric(n), dbh2 = numeric(n),
	dbh3 = numeric(n), dbh4 = numeric(n), dbh5 = numeric(n))
ls_plot_id = unique(treeData[, pointInventory_id])

clim_ind = indices[type == "parent", index_clim_start]

for (t in 1:5)
{
	precip = (climate[clim_ind, pr] - mu_pr)/sd_pr
	current_col = paste0("dbh", t - 1)
	current_dbh = latent_dbh[, ..current_col]
	setnames(current_dbh, new = "lala")
	current_dbh = current_dbh[, lala]
	latent_dbh[, (paste0("dbh", t)) := current_dbh + rgamma(n, shape = (growth_fct(current_dbh, precip, params))^2/processError,
		rate = (growth_fct(current_dbh, precip, params))/processError)]
	clim_ind = clim_ind + 1
}

## Replace the observations by the simulated dummy data; do not forget the measurement error!
treeData[parents_index, dbh := rnorm(.N, mean = latent_dbh[, dbh0], sd = measurementError)]
treeData[children_index, dbh := rnorm(.N, mean = latent_dbh[, dbh5], sd = measurementError)]

any(treeData[children_index, dbh] - treeData[parents_index, dbh] < 0)
(sd_dbh = treeData[, sd(dbh)]) # Do not forget to update it to the dummy data!

#### Stan model
## Define stan variables
# Common variables
maxIter = 1200
n_chains = 3

# Initial value for states only
average_G = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

range(average_G)

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen], normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]$latent_dbh_parents))

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

	# Explanatory variables
	precip = climate[, pr], # Annual precipitations (sum over 12 months)
	pr_mu = climate_mu_sd[variable == "pr", mu],
	pr_sd = climate_mu_sd[variable == "pr", sd],

	totalTrunkArea = treeData[, totalTrunkArea]
)

## Compile model
model = cmdstan_model("./growth_2ndOption.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 10, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 13, adapt_delta = 0.9)

results$save_output_files(dir = paste0("./", species, "/"), basename = "test_procErrorFixed_prior0.1",
	timestamp = FALSE, random = FALSE)
results$save_object(file = paste0("./", species, "/test_procErrorFixed_prior0.1.rds"))

#### Check-up
## Flags
results$cmdstan_diagnose()
results$print()

## Check dbh of individual at line row < n_indiv
row = 340

latent = results$draws(paste0("latent_dbh[",
	indices[type =="parent"][row, index_gen]:indices[type =="child"][row, index_gen],
	"]"))
latent = apply(latent, 3, mean)
real = unlist(latent_dbh[row,])
measured = c(treeData[parents_index[row], dbh], treeData[children_index[row], dbh])
x = 2000:2005

ymin = min(sd_dbh*latent, real, measured)
ymax = max(sd_dbh*latent, real, measured)

op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(2000:2005, sd_dbh*latent, pch = 19, col = "#34568B", cex = 4,
	xlab = "Years (fake)", ylab = "Diameter at breast height (scaled)", ylim = c(ymin, ymax))
points(x = c(2000, 2005), y = measured, pch = 19, col = "#FA7A35", cex = 2)
points(x, real, pch = 19, col = "#CD212A", cex = 1.5)
legend(x = "topleft", legend = c("Estimated", "Measured", "Real"), fill = c("#34568B", "#FA7A35", "#CD212A"), box.lwd = 0)
dev.off()

#### Compute residuals: compare data versus latent states with obs error
n_rep = results$metadata()$iter_sampling * results$num_chains()

dt_dharma = data.table(
	rep_latent_id = c(rep(indices[type == "parent", index_gen], each = n_rep), rep(indices[type == "child", index_gen], each = n_rep)),
	rep_dbh = c(rep(treeData[parents_index, dbh], each = n_rep), rep(treeData[children_index, dbh], each = n_rep)),
	sampled = numeric(n_rep * treeData[, .N]))

dt_dharma = data.table(
	rep_latent_id = rep(indices[type == "child", index_gen], each = n_rep),
	rep_dbh = rep(treeData[children_index, dbh], each = n_rep),
	sampled = numeric(n_rep * treeData[, .N]/2))

dt_dharma = data.table(
	rep_latent_id = rep(indices[type == "parent", index_gen], each = n_rep),
	rep_dbh = rep(treeData[parents_index, dbh], each = n_rep),
	sampled = numeric(n_rep * treeData[, .N]/2))

latent_dbh_array = results$draws("latent_dbh") # dimension: iter_sampling * n_chains * number latent states
dt_dharma[, sampled := myPredictObs(latent_dbh_array, rep_latent_id), by = rep_latent_id]

sims = matrix(data = dt_dharma[, sampled], nrow = n_rep, ncol = treeData[, .N]) # each column is for one data point
sims = t(sims) # Transpose the matrix for dharma
dim(sims)

if (ncol(sims) != n_rep)
	stop("Dimensions mismatch")

if (nrow(sims) != length(unique(dt_dharma$rep_latent_id)))
	stop("Dimensions mismatch")

dt_dharma[, diff := rep_dbh/sd_dbh - sampled]
meanDiff = dt_dharma[, mean(diff), by = rep_latent_id]
print(paste0("The mean is ", round(meanDiff[, mean(V1)], digits = 8)))
plot(meanDiff[, rep_latent_id], meanDiff[, V1], xlab = "Latent tree id", ylab = "Average difference: measured - simulated", pch = 19)

index_parent_extreme = which(meanDiff[rep_latent_id %in% indices[type == "parent", index_gen], V1] > 0.01)
latent_id_parent_extreme = meanDiff[index_parent_extreme, rep_latent_id]
tree_point_id_extreme = indices[index_gen %in% latent_id_parent_extreme, .(tree_id, pointInventory_id)]
latent_id_child_extreme = integer(tree_point_id_extreme[, .N])

for (i in 1:tree_point_id_extreme[, .N])
{
	latent_id_child_extreme[i] = indices[(type == "child") & (tree_id == tree_point_id_extreme[i, tree_id]) &
		(pointInventory_id == tree_point_id_extreme[i, pointInventory_id]), (index_gen)]
}

pdf("meanDiff.pdf", height = 8, width = 8)
plot(1:meanDiff[, .N], meanDiff[, V1], pch = 19, xlab = "Index", ylab = "Average difference: measured - simulated")
abline(v = meanDiff[, .N]/2, col = "#CD212A", lwd = 4)
points(which(meanDiff[, rep_latent_id] %in% c(latent_id_parent_extreme, latent_id_child_extreme)),
	meanDiff[rep_latent_id %in% c(latent_id_parent_extreme, latent_id_child_extreme), V1], col = "#FA7A35", pch = 19, cex = 1.25)
dev.off()

hist(meanDiff[, V1])
hist(dt_dharma[, diff])

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dt_dharma[seq(1, .N, by = n_rep), rep_dbh/sd_dbh]) # treeData[parents_index, dbh/sd_dbh]

pdf("residuals_parent_test4.pdf")
plot(forDharma)
dev.off()

#### Retrieve the parameters back
# If there is no mistake (see notebook), then the relation between some stan parameters and the parameters for non standardised dbh is:
#	potentialGrowth_stan = potentialGrowth_R - log(sd(dbh))
#	slope_dbh_stan = sd(dbh)*slope_dbh_R,
# Where *_stan are for the stan estimates, while *_R are for the real values in R
# Note that there is no transformation required for the climate parameters, only the dbh_slope and the intercept are influenced






#### 

par(mfrow = c(2,1))
lazyTrace(latent_dbh_array[, , paste0("latent_dbh[", indices[type =="parent"][row, index_gen], "]")],
	val1 = treeData[parents_index[row], dbh]/sd(treeData[, dbh]), main = paste0("tree ", row, " at t0"), scaling = sd(treeData$dbh))
lazyTrace(latent_dbh_array[, , paste0("latent_dbh[", indices[type =="child"][row, index_gen], "]")],
	val1 = treeData[children_index[row], dbh]/sd(treeData[, dbh]), main = paste0("tree ", row, " at t1"), scaling = sd(treeData$dbh))
dev.off()


plot(dataGrowth[, dbh], dataGrowth[, growth], pch = 19)
curve(growth_fct(x, params, 0, sd_dbh), col = "#178F92", lwd = 3, from = 0, to = 700, add = TRUE)

lazyTrace(results$draws("potentialGrowth"))
lazyTrace(results$draws("slope_dbh"))


## Few trace plots
lazyTrace(results$draws("potentialGrowth"), val1 = potentialGrowth)
lazyTrace(results$draws("dbh_slope"), val1 = dbh_slope)
lazyTrace(results$draws("pr_slope"), val1 = pr_slope)
lazyTrace(results$draws("pr_slope2"), val1 = pr_slope2)

lazyTrace(sd_dbh^2*results$draws("processError")) # processError is a variance, so should be multiplied by sd(dbh)^2
lazyTrace(sd_dbh*results$draws("measureError")) # measureError is a standard deviation, so should be multiplied by sd(dbh)

#### Check likelihood (for second test, see notebook)
vec_pG = seq(potentialGrowth - 0.1, potentialGrowth + 0.1, length.out = 40)
vec_dbhSlope = seq(dbh_slope - 0.008, dbh_slope + 0.012, length.out = 40)

ll = matrix(data = 0, nrow = length(vec_pG), ncol = length(vec_dbhSlope))
control = vector(mode = "list", length = ncol(ll)*nrow(ll))
names(control) = paste0("iter_", 1:length(control))

dim(ll)
r = 1
for (pG in vec_pG)
{
	c = 1
	for (dbhSlope in vec_dbhSlope)
	{
		currentParams = c(potentialGrowth = pG, dbh_slope = dbhSlope, pr_slope = pr_slope, pr_slope2 = pr_slope2)
		ll[r, c] = loglik(the_answer = the_answer/treeData[, sd(dbh)], n_indiv = n_indiv, n_hiddenState = indices[.N, index_gen],
			nbYearsPerIndiv = nbYearsPerIndiv,
			climate_index = indices[type == "parent", index_clim_start],
			normalised_precip = (climate[, pr] - climate_mu_sd[1, mu])/climate_mu_sd[1, sd],
			params = currentParams, processError = processError/treeData[, var(dbh)])
		c = c + 1
	}
	if (r %% 5 == 0)
		print(r)
	r = r + 1
}

x_sol = mean(results$draws("dbh_slope"))
y_sol = mean(results$draws("potentialGrowth"))
x_real = params["dbh_slope"]
y_real = params["potentialGrowth"]

# pdf("./Tilia_platyphyllos/test_errorsFixed_latentGiven_prSlopesNull_seed=123.pdf")
filled.contour(x = vec_dbhSlope, y = vec_pG, z = t(ll), xlab = "dbh slope", ylab = "potential growth",
	plot.axes = {
		axis(1); axis(2);
		points(x_sol, y_sol, col = "#1D1F54", pch = 19, cex = 2); text(x_sol, y_sol, pos = 4, labels = "Estimated")
		points(x_real, y_real, col = "#1D7FDF", pch = 19, cex = 2); text(x_real, y_real, pos = 4, labels = "Real solution")
	})
dev.off()
