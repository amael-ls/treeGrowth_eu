
#### Aim of prog: Analysing results (check-up residuals, plots)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool function
## Get fixed values parameters
getParams = function(model_cmdstan, params_names, type = "mean")
{
	if (!(type %in% c("mean", "median")))
		stop("Unknown type. Please choose median or mean")
	
	vals = numeric(length(params_names))
	names(vals) = params_names
	for (i in 1:length(params_names))
	{
		vals[i] = ifelse(type == "mean",
			mean(model_cmdstan$draws(params_names[i])),
			median(model_cmdstan$draws(params_names[i])))
	}
	return (vals)
}

## Get name of the last run
getLastRun = function(path, begin = "growth-", extension = ".rds", format = "ymd", run = NULL, getAll = FALSE)
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

## Bayesplot is having troubles on my mac (Arial font not always found), so I create my own traces plot
lazyTrace = function(draws, filename = NULL, run = NULL, ...)
{
	if (!is.array(draws) & !all(class(draws) %in% c("draws_array", "draws", "array")))
		stop("The class of draws should be either array, or compatible with cmdstanr (draws_array, draws, array)")
	
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
	label_ind = stri_detect(str = ls_names, regex = "label")

	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	if (any(label_ind))
		par(mar = c(5, 4, 4, 4))
	
	# Plot
	if (!is.null(filename))
	{
		pdf(paste0(filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf"))
		print(paste0("Figure saved under the name:", filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf"))
	}
	
	plot(0, pch = "", xlim = c(0, n_iter), ylim = scaling*c(min_val, max_val), axes = TRUE, bg = "transparent",
		xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""),
		ylab = ifelse(any(ylab_ind), providedArgs[["ylab"]], ""),
		main = ifelse(any(main_ind), providedArgs[["main"]], ""))

	for (chain in 1:n_chains)
	{
		if (all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(1:n_iter, scaling*draws[, chain,], type = "l", col = cols[chain])
		if (is.array(draws) & !all(class(draws) %in% c("draws_array", "draws", "array")))
			lines(1:n_iter, scaling*draws[, chain], type = "l", col = cols[chain])
	}

	if (any(val_ind))
	{
		for (val in ls_names[val_ind])
			abline(h = scaling*providedArgs[[val]], col = "#CD212A", lwd = 4)
		
		if (any(label_ind))
		{
			num_vals = stri_sub(str = ls_names[val_ind], from = stri_locate(str = ls_names[val_ind], regex = "val")[, "end"] + 1)
			for (label in ls_names[label_ind])
			{
				num_label = stri_sub(str = label, from = stri_locate(str = label, regex = "label")[, "end"] + 1)
				corresponding_val = (ls_names[val_ind])[num_vals == num_label]
				axis(4, at = scaling*providedArgs[[corresponding_val]], providedArgs[[label]], las = 1)
			}
		}
	}

	if (!is.null(filename))
		dev.off()
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

## Function to plot the prior and posterior of a parameter
lazyPosterior = function(draws, fun = dnorm, filename = NULL, run = NULL, ...)
{
	# Check-up
	if (!is.array(draws))
		stop("Draws should be an array extracted from a CmdStanMCMC object")
		
	if (!isTRUE(all.equal(fun, dnorm)) & !isTRUE(all.equal(fun, dgamma))) # isFALSE will not work here, hence !isTRUE
		stop("This function only accepts dnorm or dgamma as priors")

	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)
	nbArgs = length(providedArgs)

	# Get the argument for density if provided
	n = 512
	n_ind = (ls_names == "n")
	if (any(n_ind))
	{
		n = providedArgs[["n"]]
		print(paste0("Using n = ", n, " for the density plot"))
	}

	# Get parameters for prior
	if (isTRUE(all.equal(fun, dnorm)))
	{
		if (!all(c("mean", "sd") %in% names(providedArgs)))
			stop("You must provide mean and sd for dnorm")
		
		arg1 = providedArgs[["mean"]]
		arg2 = providedArgs[["sd"]]
	}

	if (isTRUE(all.equal(fun, dgamma)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape", "rate") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dgamma")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = temp1^2/temp2 # shape
			arg2 = temp1/temp2 # rate
		}

		if (all(c("shape", "rate") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["rate"]]
		}
	}

	# Get posterior
	density_from_draws = density(draws, n = n)
	x = density_from_draws$x
	min_x = min(x)
	max_x = max(x)

	min_x = ifelse (min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
	max_x = ifelse (max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x

	# Plot
	if (!is.null(filename))
	{
		pdf(paste0(filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf"))
		print(paste0("Figure saved under the name:", filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf"))
	}
	
	# Plot posterior
	plot(density_from_draws, xlim = c(min_x, max_x), col = "#34568B", lwd = 2, main = "Prior and posterior measurement error")
	polygon(density_from_draws, col = "#34568B22")

	# Plot prior
	curve(fun(x, arg1, arg2), add = TRUE, lwd = 2, col = "#F4C430")
	DescTools::Shade(fun(x, arg1, arg2), breaks = c(min_x, max_x), col = "#F4C43022", density = NA)

	legend(x = "topleft", legend = c("Prior", "Posterior"), fill = c("#F4C430", "#34568B"), box.lwd = 0)

	if (!is.null(filename))
		dev.off()
}

#### Read data
## Common variables
species = "Picea_abies"
path = paste0("./", species, "/")
run = 1
info_lastRun = getLastRun(path = path, run = run)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]

isPrecip_normalised = TRUE
isDBH_normalised = TRUE

if (isPrecip_normalised)
{
	print("Precip parameter must be transformed when working on the real precip scale")
	norm_clim_dt = readRDS(paste0(path, "climate_normalisation.rds"))
}

if (isDBH_normalised)
{
	print("DBH must be transformed when working on the real DBH scale")
	norm_dbh_dt = readRDS(paste0(path, "dbh_normalisation.rds"))
	sd_dbh = norm_dbh_dt[, sd]
}

## Load results
results = readRDS(paste0(path, lastRun))

#### Plot prior and posterior
## For measurement error
measurementError_array = results$draws("measureError")
lazyPosterior(draws = measurementError_array, fun = dgamma, filename = paste0(path, "measureError_posterior"), run = run,
	shape = 3.0/0.001, rate = sd_dbh*sqrt(3.0)/0.001)

## For process error
processError_array = results$draws("processError")
lazyPosterior(draws = processError_array, fun = dgamma, filename = paste0(path, "processError_posterior"),
	shape = 5.0^2/1, rate = sd_dbh^2*5.0/1)

#### Plot chains main parameters
lazyTrace(draws = results$draws("potentialGrowth"), filename = paste0(path, "potentialGrowth"), run = run)
lazyTrace(draws = results$draws("dbh_slope"), filename = paste0(path, "dbh_slope"), run = run)

lazyTrace(draws = results$draws("pr_slope"), filename = paste0(path, "pr_slope"), run = run)
lazyTrace(draws = results$draws("pr_slope2"), filename = paste0(path, "pr_slope2"), run = run)
lazyTrace(draws = results$draws("tas_slope"), filename = paste0(path, "tas_slope"), run = run)
lazyTrace(draws = results$draws("tas_slope2"), filename = paste0(path, "tas_slope2"), run = run)

lazyTrace(draws = results$draws("competition_slope"), filename = paste0(path, "competition_slope"), run = run)

#### Get parameters and convert them to the 'real' scale
# If there is no mistake (see notebook), then the relation between some stan parameters and the parameters for non standardised dbh is:
#	potentialGrowth_s = potentialGrowth_r - log(sd(dbh))
#	dbh_slope_s = sd(dbh)*dbh_slope_r,
# Where *_s are for the estimates on the standardised scale, while *_r are for the estimates on the 'real' scale (i.e. non-standardised)
# Note that there is no transformation required for the climate parameters, only the dbh_slope and the intercept are influenced

(potentialGrowth = mean(results$draws("potentialGrowth")) + log(sd_dbh))
(dbh_slope = mean(results$draws("dbh_slope"))/sd_dbh)
(pr_slope = mean(results$draws("pr_slope")))
(pr_slope2 = mean(results$draws("pr_slope2")))

#### Posterior predictive checking: Can the model give rise to new observations that properly resemble the original data?
## Script to generate new data. Note that 'model' is an unnecessary block here as restart from results. gq stands for generated quantities
# More information can be found at https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html
# The new observations are simulated as follow:
#	1. Draw the vector of parameters theta (which includes the latent states!)
#	2. Generate the parent observation from the corresponding latent dbh according to the model
#	3. Generate the children observation according to the model (i.e., after n years of yearly latent growth)

gq_script = write_stan_file(
"
functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real totalTreeWeight,
		real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0>[n_indiv] totalTreeWeight; // Sum of the tree weights for a given plot at a given tieme
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_indiv] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);
	// Normalised and centred totalTreeWeight
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd(Yobs)^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	// Declaration output variables
	array [n_obs] real newObservations;
	array [n_obs] real latent_dbh_parentsChildren;
	array [n_latentGrowth] real procError_sim;
	array [n_latentGrowth + n_indiv] real yearly_latent_dbh;
	
	{
		// Variables declared in nested blocks are local variables, not generated quantities, and thus won't be printed.
		real current_latent_dbh;
		int growth_counter = 1;
		int dbh_counter = 1;
		real expected_growth;

		for (i in 1:n_indiv) // Loop over all the individuals
		{
			// Fill the parents
			latent_dbh_parentsChildren[parents_index[i]] = latent_dbh_parents[i];
			yearly_latent_dbh[dbh_counter] = latent_dbh_parents[i];

			// Generate the parent observation conditional on the parent state
			newObservations[parents_index[i]] = normal_rng(latent_dbh_parents[i], measureError);

			// Starting point to compute latent dbh child
			current_latent_dbh = latent_dbh_parents[i];

			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1], normalised_totalTreeWeight[i],
					potentialGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);

				// Record difference expected growth versus latent growth
				procError_sim[growth_counter] = expected_growth - latent_growth[growth_counter];

				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter]; // Or should it be += gamma(mean = latent_growth, var = processError)?
				growth_counter += 1;
				dbh_counter += 1;
				yearly_latent_dbh[dbh_counter] = current_latent_dbh;
			}

			// The last current_latent_dbh corresponds to the child latent dbh
			latent_dbh_parentsChildren[children_index[i]] = current_latent_dbh;
			newObservations[children_index[i]] = normal_rng(latent_dbh_parentsChildren[children_index[i]], measureError);
			dbh_counter += 1;
		}
	}
}
"
)

## Compile simulation generator
gq_model = cmdstan_model(gq_script)

## Generate simulations
# Access data
n_chains = results$num_chains()
iter_sampling = results$metadata()$iter_sampling
stanData = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
n_obs = stanData$n_obs
n_hiddenState = stanData$n_latentGrowth + stanData$n_indiv
indices = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "indices.rds"))

# Simulations
generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)
dim(generate_quantities$draws()) # iter_sampling * n_chains * (2*n_obs + n_latentGrowth + n_hiddenState)

## Check that the observation residuals are ok. There should not be any difference between parents and children, and no pattern with dbh
n_rep = iter_sampling * n_chains

dt_dharma = data.table(rep_id = rep(1:n_obs, each = n_rep),
	rep_dbh = rep(stanData$Yobs, each = n_rep), # Observations
	latent_dbh = numeric(n_rep * n_obs), # Estimated latent dbh
	simulated_observations = numeric(n_rep * n_obs))

newObservations_array = generate_quantities$draws("newObservations")
latent_dbh_array = generate_quantities$draws("latent_dbh_parentsChildren")

dt_dharma[, simulated_observations := myPredictObs(newObservations_array, rep_id, regex = "newObservations"), by = rep_id]
dt_dharma[, simulated_observations := sd_dbh*simulated_observations]

dt_dharma[, latent_dbh := myPredictObs(latent_dbh_array, rep_id, regex = "latent_dbh_parentsChildren"), by = rep_id]
dt_dharma[, latent_dbh := sd_dbh*latent_dbh]

dt_dharma[, residuals_obs := rep_dbh - simulated_observations]

dt_dharma[rep_id %in% stanData$parents_index, type := "parents"]
dt_dharma[rep_id %in% stanData$children_index, type := "children"]

setorderv(x = dt_dharma, cols = c("type", "rep_id"), order = c(-1, 1)) # The -1 is to have parents first

jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_scatter.jpg"), quality = 50)
plot(dt_dharma[, rep_id], dt_dharma[, residuals_obs], pch = 19, col = "#A1A1A122")
abline(v = stanData$n_indiv, lwd = 2, col = "#CD212A")
axis(3, at = n_obs/4, "Parents", las = 1)
axis(3, at = 3*n_obs/4, "Children", las = 1)
dev.off()

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_hist.pdf"))
hist(dt_dharma[, residuals_obs])
dev.off()

qq = qqnorm(dt_dharma[, residuals_obs], pch = 1, frame = FALSE, plot.it = FALSE)

jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "qqplot.jpg"))
plot(qq)
qqline(dt_dharma[, residuals_obs], col = "#34568B", lwd = 3)
dev.off()

## Check that the process error does not show any pattern with any predictor, and that it is gamma distributed
processError_array = generate_quantities$draws("procError_sim")
mean(processError_array)

dim(processError_array) # iter_sampling * n_chains * n_latentGrowth

processError_vec = as.vector(processError_array) # The order is array[,1,1], array[,2,1], ..., array[,n_chain,1], array[,2,1], ...
length(processError_vec)

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "processError_check_hist.pdf"))
hist(processError_vec)
dev.off()

processError_avg = apply(processError_array, 3, mean) # Average processError for each latent growth (or dbh)
length(processError_avg)

index_notLastMeasure = 1:n_hiddenState
index_notLastMeasure = index_notLastMeasure[!(index_notLastMeasure %in% indices[stanData$children_index, index_gen])]

if (length(index_notLastMeasure) != length (processError_avg))
	stop("Dimension mismatch between processError_avg and index_notLastMeasure")

latent_dbh = apply(generate_quantities$draws("yearly_latent_dbh"), 3, mean)
latent_dbh = latent_dbh[index_notLastMeasure]

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "processError_vs_dbh_check.pdf"))
plot(latent_dbh, processError_avg, type = "l")
dev.off()

cor(latent_dbh, processError_avg)











####! CRASH TEST ZONE
#! TEST
ll = bayesplot::nuts_params(results)
test = bayesplot::mcmc_nuts_energy(ll)
marginalPlot(results$draws("measureError"))
lazyTrace(results$draws("measureError", inc_warmup = TRUE))
setDT(ll)
ll = dcast(ll, Chain + Iteration ~ Parameter, value.var = "Value")
#! END TEST



#### Get parameters
paramsNames = c("potentialMaxGrowth", "power_dbh", "optimal_precip", "width_precip_niche", "processError")

## Mean values
meanParams = getParams(results, paramsNames)
processError_mean = meanParams[["processError"]]

## Median values
medParams = getParams(results, paramsNames, type = "median")
processError_med = medParams[["processError"]]

#### Residuals
## Load data
treeFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
# treeData = readRDS(paste0(treeFolder, "trees_forest_reshaped.rds"))
treeData = readRDS(paste0(path, "trees_forest_reshaped.rds")) # Folder to use when running on local computer
treeData = treeData[speciesName_sci == species]

## Get dbh and time
dbh_start = treeData[treeData[, .I[which.min(year)], by = .(tree_id, pointInventory_id)][, V1], dbh]
dbh_end = treeData[treeData[, .I[which.max(year)], by = .(tree_id, pointInventory_id)][, V1], dbh]

t_start = treeData[treeData[, .I[which.min(year)], by = .(tree_id, pointInventory_id)][, V1], year]
t_end = treeData[treeData[, .I[which.max(year)], by = .(tree_id, pointInventory_id)][, V1], year]
delta_t = t_end - t_start

## Climate
climFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/climateData/Chelsa/yearlyAverage/"
# climate = readRDS(paste0(climFolder, "FR_reshaped_climate.rds"))
climate = readRDS(paste0(path, "FR_reshaped_climate.rds")) # Folder to use when running on local computer

## indices
# indices = readRDS(paste0(treeFolder, species, "_indices.rds"))
indices = readRDS(paste0(path, species, "_indices.rds"))
indices_data = indices[, index_gen] 
indices = unique(indices[, .(tree_id, pointInventory_id, index_clim_start, index_clim_end)])

## Create simulations
n_rep = 250
dt = data.table(rep_dbh_end = rep(dbh_end, each = n_rep), sampled = numeric(n_rep * length(dbh_end)))
medPred = numeric(length(dbh_end))

if (isDBH_normalised)
	dt[, rep_dbh_end := rep_dbh_end/norm_dbh_dt[, sd]]

i = 1
yr = 1
for (i in 1:length(dbh_end))
{
	currentDbh = dbh_start[i]/ifelse(isDBH_normalised, norm_dbh_dt[, sd], 1)
	currentDbh_med = dbh_start[i]/ifelse(isDBH_normalised, norm_dbh_dt[, sd], 1)
	currentClim_index = indices[i, index_clim_start]
	for (yr in 1:delta_t[i])
	{
		precip = climate[currentClim_index, pr]
		currentDbh = rnorm(n_rep, rnorm(n_rep, nextDBH(currentDbh, precip, meanParams, norm_clim_dt), processError_mean), 0.006)
		currentDbh_med = nextDBH(currentDbh_med, precip, medParams, norm_clim_dt)
		currentClim_index = currentClim_index + 1
	}
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := currentDbh]
	medPred[i] = currentDbh_med
	if (i %% 100 == 0)
		print(paste0(round(i*100/length(dbh_end), 2), "% done"))
}

# saveRDS(dt, "testForResiduals.rds")

## Reshape simulations into a matrix length(dbh) x n_rep
sims = matrix(data = dt[, sampled], nrow = n_rep, ncol = length(dbh_end)) # each column is for one data point (the way I organised the dt)
sims = t(sims) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dbh_end/ifelse(isDBH_normalised, norm_dbh_dt[, sd], 1)
, fittedPredictedResponse = medPred)

plot(forDharma)
dev.off()

# #### Residuals using posterior samples directly (should have nIter/2 samples per measure)
# ## Get the posterior distrib of the states that have data to be compared to
# n_rep = results$metadata()$iter_sampling
# dt = data.table(measured = rep(treeData[, dbh], each = n_rep), sampled_mean = numeric(n_rep * treeData[, .N]),
# 	sampled_med = numeric(n_rep * treeData[, .N]))
# if (isDBH_normalised)
# 	dt[, measured := measured/norm_dbh_dt[, sd]]

# for (i in 1:treeData[, .N])
# {
# 	dt[((i - 1)*n_rep + 1):(i*n_rep),
# 		sampled_mean := rnorm(n_rep, rowMeans(results$draws(paste0("latent_dbh[", indices_data[i], "]"))), 0.006)]
# 	dt[((i - 1)*n_rep + 1):(i*n_rep),
# 		sampled_med := rnorm(n_rep, apply(results$draws(paste0("latent_dbh[", indices_data[i], "]")), 1, median), 0.006)]
# 	if (i %% 200 == 0)
# 		print(paste0(round(i*100/treeData[, .N], 2), "% done"))
# }

# ## Reshape simulations into a matrix length(dbh) x n_rep
# sims = matrix(data = dt[, sampled_mean], nrow = n_rep, ncol = treeData[, .N]) # each column is for one data point (the way I organised the dt)
# sims = t(sims) # Transpose the matrix for dharma

# forDharma = createDHARMa(simulatedResponse = sims,
# 	observedResponse = dt[seq(1, .N, by = n_rep), measured]) # all.equal(ll, treeData[, dbh]/norm_dbh_dt[, sd]) 
# # , fittedPredictedResponse = dt[seq(1, .N, by = n_rep), sampled_med])

# # fittedPredictedResponse: optional fitted predicted response. For
# #	Bayesian posterior predictive simulations, using the median
# #	posterior prediction as fittedPredictedResponse is
# #	recommended. If not provided, the mean simulatedResponse will
# #	be used.

# plot(forDharma)
# dev.off()

# plot(dt[, measured], dt[, sampled])

#### Plot posterior distributions and traces of parameters
# Folder to save plots
figurePath = paste0(path, time_ended, "/")
if (!dir.exists(figurePath))
	dir.create(figurePath)

# Potential maximal growth
plot_title = ggplot2::ggtitle("Posterior distribution potentialMaxGrowth", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("potentialMaxGrowth"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "pmg-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for potentialMaxGrowth")
temp_plot = mcmc_trace(results$draws("potentialMaxGrowth")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "pmg-traces.pdf"), plot = temp_plot, device = "pdf")

# Power dbh
plot_title = ggplot2::ggtitle("Posterior distribution power_dbh", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("power_dbh"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "power_dbh-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for power_dbh")
temp_plot = mcmc_trace(results$draws("power_dbh")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "power_dbh-traces.pdf"), plot = temp_plot, device = "pdf")

# Optimal precipitation
plot_title = ggplot2::ggtitle("Posterior distribution optimal_precip", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("optimal_precip"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "optimal_precip-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for optimal_precip")
temp_plot = mcmc_trace(results$draws("optimal_precip")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "optimal_precip-traces.pdf"), plot = temp_plot, device = "pdf")

# Width niche distribution (precipitation)
plot_title = ggplot2::ggtitle("Posterior distribution width_precip_niche", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("width_precip_niche"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "width_precip_niche-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for width_precip_niche")
temp_plot = mcmc_trace(results$draws("width_precip_niche")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "width_precip_niche-traces.pdf"), plot = temp_plot, device = "pdf")

# processError
plot_title = ggplot2::ggtitle("Posterior distribution processError", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("processError"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "processError-distrib.pdf"), plot = temp_plot, device = "pdf")

lazyTrace(results$draws("processError"))

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

lazyTrace(results$draws("measureError"))

cor(results$draws("processError"), results$draws("measureError"))

for (i in c(1:6, 12:14))
{
	plot_title = ggplot2::ggtitle(paste0("Posterior distribution latent_dbh[", i, "]"),
		"with medians and 80% intervals")
	temp_plot = mcmc_areas(results$draws(paste0("latent_dbh[", i, "]")), prob = 0.8) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latent_dbh_", i, "-distrib.pdf"), plot = temp_plot, device = "pdf")

	plot_title = ggplot2::ggtitle(paste0("Traces for latent_dbh[", i, "]"))
	temp_plot = mcmc_trace(results$draws(paste0("latent_dbh[", i, "]"), inc_warmup = TRUE)) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latent_dbh_", i, "-traces.pdf"), plot = temp_plot, device = "pdf")
}


latent_1_6 = getParams(results, paste0("latent_dbh[", 1:6, "]"))
x = 2000:2005

pdf("./latent_real_2ndOption.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(2000:2005, norm_dbh_dt[, sd]*latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)")
points(x = 2000, y = treeData[1, dbh], pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = treeData[2, dbh], pch = 19, col = "#CD212A", cex = 2)
dev.off()

pdf("./growth_dbh.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
curve(norm_dbh_dt[, sd]^(1 - meanParams["power_dbh"]) * growth(x, 0, params = meanParams, norm_clim_dt, isScaled = TRUE),
	from = 50, to = 1000, lwd = 2, col = "#34568B", xlab = "dbh", ylab = "Growth (in mm/yr)")
dev.off()

# pdf("./growth_precip.pdf", height = 7, width = 7)
# op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
# curve(growth(250/norm_dbh_dt[, sd], x, params = meanParams, norm_clim_dt, isScaled = FALSE), from = 650, to = 1500,
# 	lwd = 2, col = "#34568B", xlab = "Precipitation (mm/yr)", ylab = "Growth (in mm/yr)")
# dev.off()

## Compute average climate for each tree during its time span
indices[, avg_pr := mean(climate[index_clim_start:index_clim_end, pr]), by = .(tree_id, pointInventory_id)]

# What follows works only for the French data that have a particular structure. For the general case, use indices
treeData = treeData[indices[, .(tree_id, pointInventory_id, avg_pr)], on = c("tree_id", "pointInventory_id")]
dbh0 = treeData[seq(1, .N - 1, by = 2), dbh]
dbh1 = treeData[seq(2, .N, by = 2), dbh]
yearly_avg_growth = (dbh1 - dbh0)/(treeData[seq(2, .N, by = 2), year] - treeData[seq(1, .N - 1, by = 2), year])

mean(treeData[, avg_pr])
split_clim = min(treeData[, avg_pr]) + 0:3*(max(treeData[, avg_pr]) - min(treeData[, avg_pr]))/3

treeData[, class_pr := ""]
treeData[(split_clim[1] <= avg_pr) & (avg_pr < split_clim[2]), class_pr := "low"]
treeData[(split_clim[2] <= avg_pr) & (avg_pr < split_clim[3]), class_pr := "medium"]
treeData[(split_clim[3] <= avg_pr) & (avg_pr <= split_clim[4]), class_pr := "high"]

colours = c("#CD212A", "#DAAE29", "#094F9A")
treeData[class_pr == "low", colour := colours[1]]
treeData[class_pr == "medium", colour := colours[2]]
treeData[class_pr == "high", colour := colours[3]]

# jpeg("test.jpg", width = 1080, height = 1080, quality = 100)
pdf("test.pdf", width = 7, height = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(dbh0, yearly_avg_growth, pch = 19, col = treeData[seq(1, .N - 1, by = 2), colour],
	xlab = "dbh (at time 0)", ylab = "yearly average growth", cex = 0.35)
dev.off()

####! END CRASH TEST ZONE

