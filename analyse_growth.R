
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
	{
		print(paste("Searching among runs =", run))
		begin = paste0(begin, "run=", run, "-")
	}
	
	ls_files = list.files(path = path, pattern = paste0("^", begin, ".*", extension, "$"))

	if (length(ls_files) == 0)
	{
		warning(paste0("No file detected in the folder '", path, "'. You were looking for '", begin, "*", extension, "'"))
		return (list(file = NA, time_ended = NA))
	}

	if (is.null(run))
		ls_files = ls_files[!stri_detect(str = ls_files, regex = paste0(begin, "run="))]

	ls_files_split = stri_split(
		str = stri_sub(str = ls_files,
			from = stri_locate(ls_files, regex = begin)[, "end"] + 1,
			to = stri_locate_last(ls_files, regex = "_[[:digit:]].*.rds")[, "start"] - 1),
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

## Function to reshape draws_array
reshapeDraws = function(draws_array, id_latent, regex = "latent_dbh")
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
		
	# isFALSE will not work here, hence !isTRUE
	if (!isTRUE(all.equal(fun, dnorm)) & !isTRUE(all.equal(fun, dgamma)) & !isTRUE(all.equal(fun, dbeta)))
		stop("This function only accepts dnorm, dgamma, or dbeta as priors")

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

	params = ""
	params_ind = (ls_names == "params")
	if (any(params_ind))
		params = providedArgs[["params"]]

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

	if (isTRUE(all.equal(fun, dbeta)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape1", "shape2") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape1 and shape2 for dbeta")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = ((1 - temp1)/temp2 - 1/temp1)*temp1 ^ 2 # shape 1
			arg2 = arg1*(1/temp1 - 1) # shape 2
		}

		if (all(c("shape1", "shape2") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape1"]]
			arg2 = providedArgs[["shape2"]]
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
	plot(density_from_draws, xlim = c(min_x, max_x), col = "#34568B", lwd = 2, main = paste("Prior and posterior", params))
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
species = "Fagus sylvatica" # "Tilia platyphyllos"
path = paste0("./", species, "/")
run = 3
info_lastRun = getLastRun(path = path, run = run)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]

# isClimate_normalised = TRUE
isDBH_normalised = TRUE

# if (isClimate_normalised)
# {
# 	print("Climate parameter must be transformed when working on the real climate scale")
# 	norm_clim_dt = readRDS(paste0(path, ifelse(is.null(run), "", paste0(run, "_")), "climate_normalisation.rds"))
# }

if (isDBH_normalised)
{
	print("DBH must be transformed when working on the real DBH scale")
	norm_dbh_dt = readRDS(paste0(path, ifelse(is.null(run), "", paste0(run, "_")), "dbh_normalisation.rds"))
	sd_dbh = norm_dbh_dt[, sd]
}

## Load results and associated data set
results = readRDS(paste0(path, lastRun))
stanData = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))

#### Plot prior and posterior
## For process error (sigmaProc), it is a variance (of a gamma distrib)
sigmaProc_array = results$draws("sigmaProc")
lazyPosterior(draws = sigmaProc_array, fun = dgamma, filename = paste0(path, "sigmaProc_posterior"), params = "process error",
	shape = 5.0^2/1, rate = sd_dbh^2*5.0/1, run = run)

print(paste("The sd of growth is +/-", round(sd_dbh*sqrt(mean(sigmaProc_array)), 3), "mm"))

## For measurement error
for (i in 1:stanData$n_inventories)
{
	# For routine measurement error (sigmaObs), it is a sd (of a normal distrib)
	sigmaObs_array = results$draws(paste0("sigmaObs[", i, "]"))
	lazyPosterior(draws = sigmaObs_array, fun = dgamma, filename = paste0(path, "sigmaObs_posterior[", i, "]"), run = run,
		shape = 3.0/0.15, rate = sd_dbh*sqrt(3.0)/0.15, params = paste0("routine obs error[", i, "]"))

	print(paste("The routine error is", round(sd_dbh*mean(sigmaObs_array), 3), "mm."))

	# Extreme error (etaObs)...
	etaObs_array = results$draws(paste0("etaObs[", i, "]"))
	lazyPosterior(draws = etaObs_array, fun = dgamma, filename = paste0(path, "etaObs_posterior[", i, "]"), run = run,
		params = paste0("extreme obs error[", i, "]"), shape = 25.6^2/6.2, rate = sd_dbh*25.6/6.2)

	# ... and its associated probability of occurrence
	proba_array = results$draws(paste0("p[", i, "]")) #! Renamed proba, to change in the future !!!
	lazyPosterior(draws = proba_array, fun = dbeta, filename = paste0(path, "proba_posterior[", i, "]"), run = run,
		params = paste0("probability extreme obs error[", i, "]"), shape1 = 48.67, shape2 = 1714.84)

	print(paste("The extreme error is", round(sd_dbh*mean(etaObs_array), 3), "mm. It occurs with a probability p =",
		round(mean(proba_array), 4)))

	## Correlation routine obs and extreme obs
	print(paste("Correlation routine <---> extreme = ", round(cor(sigmaObs_array, etaObs_array), 3)))
}

#### Plot chains main parameters
lazyTrace(draws = results$draws("averageGrowth_mu", inc_warmup = FALSE), filename = paste0(path, "averageGrowth_mu"), run = run)
lazyTrace(draws = results$draws("averageGrowth_sd", inc_warmup = FALSE), filename = paste0(path, "averageGrowth_sd"), run = run)
lazyTrace(draws = results$draws("dbh_slope", inc_warmup = FALSE), filename = paste0(path, "dbh_slope"), run = run)

lazyTrace(draws = results$draws("pr_slope", inc_warmup = FALSE), filename = paste0(path, "pr_slope"), run = run)
lazyTrace(draws = results$draws("pr_slope2", inc_warmup = FALSE), filename = paste0(path, "pr_slope2"), run = run)
lazyTrace(draws = results$draws("tas_slope", inc_warmup = FALSE), filename = paste0(path, "tas_slope"), run = run)
lazyTrace(draws = results$draws("tas_slope2", inc_warmup = FALSE), filename = paste0(path, "tas_slope2"), run = run)

lazyTrace(draws = results$draws("competition_slope", inc_warmup = FALSE), filename = paste0(path, "competition_slope"), run = run)

# results$print(c("lp__", "averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	# "ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"))

results$print(c("lp__", "averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"competition_slope", "sigmaObs", "etaObs", "p", "sigmaProc"), max_rows = 20)

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
	real growth(real dbh0, real precip, real temperature /*, real ph*/, real standBasalArea, real averageGrowth, real dbh_slope,
		real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2 /* , real ph_slope, real ph_slope2 */, real competition_slope)
	{
		return (exp(averageGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 /* + ph_slope*ph + ph_slope2*ph^2 */ + competition_slope*standBasalArea));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Total number of individuals (all NFIs together)
	int<lower = 1> n_climate; // Dimension of the climate vector (all NFIs together)
	int<lower = 1, upper = n_indiv> n_plots; // Number of plots (all NFIs together)
	int<lower = 1> n_obs; // Total number of tree observations (all NFIs together)
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth (all NFIs together)
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children tree observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the 'observation space'
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the 'observation space'
	array [n_children] int<lower = 2> latent_children_index; // Index of children in the 'latent space'
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent
	
	array [n_indiv] int<lower = 1, upper = n_inventories> nfi_id; // Number of the NFI for a given individual
	
	array [n_indiv] int<lower = 1, upper = n_plots> plot_index; // Indicates to which plot individuals belong to

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real<lower = 0> pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	// vector<lower = 0, upper = 14>[n_plots] ph; // pH of the soil measured with CaCl2
	// real<lower = 0, upper = 14> ph_mu; // To standardise the pH
	// real<lower = 0> ph_sd; // To standardise the pH

	vector<lower = 0>[n_indiv] standBasalArea; // Sum of the tree basal area for a given plot at a given time
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	// vector[n_plots] normalised_ph = (ph - ph_mu)/ph_sd; // Normalised and centred pH
	vector[n_indiv] normalised_standBasalArea = (standBasalArea - mean(standBasalArea))/sd(standBasalArea); // Normalised and centred BA
}

parameters {
	// Parameters for growth function
	array [n_plots] real averageGrowth; // Growth (grouped by plots) when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	// real ph_slope;
	// real ph_slope2;

	real competition_slope;

	// Hyper parameters for growth function
	real averageGrowth_mu;
	real<lower = 0> averageGrowth_sd;

	// Errors (observation and process)
	// --- Process error, which is the variance of a gamma distrib /!\
	real<lower = 0.5/sd_dbh^2> sigmaProc;

	// --- Routine observation error, which is constrained by default, see appendix D Eitzel for the calculus.
	array [n_inventories] real<lower = 0.1/sqrt(12)*25.4/sd_dbh> sigmaObs; // Std. Dev. of a normal distrib /!\

	// --- Extreme error, by default at least twice the observation error. Rüger 2011 found it is around 8 times larger
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> p; // probabilities of occurrence of extreme errors (etaObs) for each NFI

	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
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
		int children_counter = 1;
		real expected_growth;

		for (i in 1:n_indiv) // Loop over all the individuals
		{
			// Fill the parents
			latent_dbh_parentsChildren[parents_index[i]] = latent_dbh_parents[i];
			yearly_latent_dbh[dbh_counter] = latent_dbh_parents[i];

			// Generate the parent observation conditional on the parent state
			newObservations[parents_index[i]] = (1 - p[nfi_id[i]]) * normal_rng(latent_dbh_parents[i], sigmaObs[nfi_id[i]]) +
				p[nfi_id[i]] * normal_rng(latent_dbh_parents[i], etaObs[nfi_id[i]]);

			// Starting point to compute latent dbh child
			current_latent_dbh = latent_dbh_parents[i];

			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1] /*, normalised_ph[plot_index[i]]*/, normalised_standBasalArea[i],
					averageGrowth[plot_index[i]], dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2 /*, ph_slope, ph_slope2 */,
					competition_slope);

				// Record difference expected growth versus latent growth
				procError_sim[growth_counter] = expected_growth - latent_growth[growth_counter];

				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter]; // Or should it be += gamma(mean = latent_growth, var = sigmaProc)?
				growth_counter += 1;
				dbh_counter += 1;
				yearly_latent_dbh[dbh_counter] = current_latent_dbh;

				// Recording the children dbh
				if (dbh_counter == latent_children_index[children_counter])
				{
					latent_dbh_parentsChildren[children_index[i]] = current_latent_dbh;
					newObservations[children_index[i]] = (1 - proba[nfi_id[i]]) * normal_rng(current_latent_dbh, sigmaObs[nfi_id[i]]) +
						proba[nfi_id[i]] * normal_rng(current_latent_dbh, etaObs[nfi_id[i]]);
					children_counter += 1;
				}
			}
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
n_obs = stanData$n_obs
n_hiddenState = stanData$n_latentGrowth + stanData$n_indiv
indices = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "indices.rds"))

indices[, nfi_index := 1]
indices[stri_detect_regex(plot_id, "germany"), nfi_index := 2]

indices[stri_detect_regex(plot_id, "france"), .N]
indices[stri_detect_regex(plot_id, "germany"), .N]

growth_counter = 1
dbh_counter = 1
children_counter = 1

latent_children_index = indices[type == "child", index_gen]
indices[, type := "child"]
indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]
length(latent_children_index)

for (i in 1:n_indiv)
{
	for (j in 1:nbYearsGrowth[i])
	{
		growth_counter = growth_counter + 1;
		dbh_counter = dbh_counter + 1;

		if (dbh_counter == latent_children_index[children_counter])
			children_counter = children_counter + 1;
	}
	dbh_counter = dbh_counter + 1;
}

growth_counter - stanData$n_latentGrowth
dbh_counter - stanData$n_latentGrowth - n_indiv
children_counter - length(latent_children_index)


stanData$nfi_id = unique(indices[, .(plot_id, tree_id, nfi_index)])[, nfi_index]

if (length(stanData$nfi_id) != stanData$n_indiv)
	stop("Dimensions mismatch")

# Simulations
generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)
dim(generate_quantities$draws()) # iter_sampling * n_chains * (2*n_obs + n_latentGrowth + n_hiddenState)
# 1000 x 3 x 172 835

## Check that the observation residuals are ok. There should not be any difference between parents and children, and no pattern with dbh
n_rep = iter_sampling * n_chains

dt_dharma = data.table(rep_id = rep(1:n_obs, each = n_rep),
	rep_dbh = rep(stanData$Yobs, each = n_rep), # Observations
	latent_dbh = numeric(n_rep * n_obs), # Estimated latent dbh
	simulated_observations = numeric(n_rep * n_obs))

newObservations_array = generate_quantities$draws("newObservations")
dim(newObservations_array) # iter_sampling * n_chains * n_obs

chain1 = newObservations_array[, 1, ] # iter_sampling * n_obs 
obs_1_chain_1 = newObservations_array[, 1, 14885] # iter_sampling
sum(is.na(obs_1_chain_1))

latent_dbh_array = generate_quantities$draws("latent_dbh_parentsChildren")
dim(latent_dbh_array) # iter_sampling * n_chains * n_obs

# There area NAs from 44 652 001 to 52 971 000, length = 4 161 000

dt_dharma[, simulated_observations := reshapeDraws(newObservations_array, rep_id, regex = "newObservations"), by = rep_id]
dt_dharma[, simulated_observations := sd_dbh*simulated_observations]

dt_dharma[, latent_dbh := reshapeDraws(latent_dbh_array, rep_id, regex = "latent_dbh_parentsChildren"), by = rep_id]
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
sigmaProc_array = generate_quantities$draws("procError_sim")
mean(sigmaProc_array)

dim(sigmaProc_array) # iter_sampling * n_chains * n_latentGrowth

sigmaProc_vec = as.vector(sigmaProc_array) # The order is array[,1,1], array[,2,1], ..., array[,n_chain,1], array[,2,1], ...
length(sigmaProc_vec)

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "sigmaProc_check_hist.pdf"))
hist(sigmaProc_vec)
dev.off()

sigmaProc_avg = apply(sigmaProc_array, 3, mean) # Average sigmaProc for each latent growth (or dbh)
length(sigmaProc_avg)

index_notLastMeasure = 1:n_hiddenState
index_notLastMeasure = index_notLastMeasure[!(index_notLastMeasure %in% indices[stanData$children_index, index_gen])]

if (length(index_notLastMeasure) != length (sigmaProc_avg))
	stop("Dimension mismatch between sigmaProc_avg and index_notLastMeasure")

latent_dbh = apply(generate_quantities$draws("yearly_latent_dbh"), 3, mean)
latent_dbh = latent_dbh[index_notLastMeasure]

jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "sigmaProc_vs_dbh_check.jpg"), quality = 50)
plot(latent_dbh, sigmaProc_avg, type = "l")
dev.off()

cor(latent_dbh, sigmaProc_avg)



#############################


gq_script = write_stan_file(
"
functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real totalTreeWeight,
		real averageGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		return (exp(averageGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
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
	real averageGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd(Yobs)^2> sigmaProc;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> sigmaObs; // Constrained by default, see appendix D Eitzel for the calculus

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
			newObservations[parents_index[i]] = normal_rng(latent_dbh_parents[i], sigmaObs);

			// Starting point to compute latent dbh child
			current_latent_dbh = latent_dbh_parents[i];

			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1], normalised_totalTreeWeight[i],
					averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);

				// Record difference expected growth versus latent growth
				procError_sim[growth_counter] = expected_growth - latent_growth[growth_counter];

				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter]; // Or should it be += gamma(mean = latent_growth, var = sigmaProc)?
				growth_counter += 1;
				dbh_counter += 1;
				yearly_latent_dbh[dbh_counter] = current_latent_dbh;
			}

			// The last current_latent_dbh corresponds to the child latent dbh
			latent_dbh_parentsChildren[children_index[i]] = current_latent_dbh;
			newObservations[children_index[i]] = normal_rng(latent_dbh_parentsChildren[children_index[i]], sigmaObs);
			dbh_counter += 1;
		}
	}
}
"
)

