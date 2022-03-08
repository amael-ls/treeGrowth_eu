
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
	requiredArgs = c("dbh_parents", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh_parents and normalise")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	normalise = providedArgs[["normalise"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	n_indiv = length(dbh_parents)
	Y_gen = rgamma(n_indiv, dbh_parents^2, dbh_parents)
	
	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
	}
	return(list(latent_dbh_parents = Y_gen))
}

## Function to rebuild latent dbh from latent dbh parents and latent growth
compute_latent_dbh = function(latent_growth_array, latent_dbh_parents_array, indices, n_hiddenState, regex = "latent_dbh",
	iter_sampling = results$metadata()$iter_sampling, num_chains = results$num_chains())
{
	n_indiv = dim(latent_dbh_parents_array)[3]
	latent_dbh_array = array(data = NA, dim = c(iter_sampling, num_chains, n_hiddenState),
		dimnames = list(1:iter_sampling, paste0("chain", 1:3), paste0(regex, "[", 1:n_hiddenState, "]")))
	parents = indices[type == "parent", index_gen]
	
	not_parent_index = 1:indices[.N, index_gen]
	not_parent_index = not_parent_index[!(not_parent_index %in% indices[type == "parent", index_gen])]

	for (i in 1:n_indiv)
		latent_dbh_array[, , paste0("latent_dbh[", parents[i], "]")] = latent_dbh_parents_array[, , paste0("latent_dbh_parents[", i, "]")]
	
	for (i in 1:(n_hiddenState - n_indiv))
	{
		ind = not_parent_index[i]
		temporary = array(latent_growth_array[, , paste0("latent_growth[", i, "]")], dim = c(iter_sampling, num_chains))
		latent_dbh_array[, , paste0("latent_dbh[", ind, "]")] = latent_dbh_array[, , paste0("latent_dbh[", ind - 1, "]")] +
			temporary
	}
	
	return (latent_dbh_array)
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

# Compute the number of growing years per individual
indices[, nbYearsGrowth := max(year) - min(year), by = .(pointInventory_id, tree_id)]

checkUp = all(indices[, nbYearsGrowth == index_clim_end - index_clim_start])
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

nbYearsGrowth = unique(indices[, .(tree_id, pointInventory_id, nbYearsGrowth)])[, nbYearsGrowth]
if (length(nbYearsGrowth) != n_indiv)
	stop("Dimension mismatch between nbYearsGrowth and n_indiv")

n_hiddenState = indices[.N, index_gen]

parents_index = treeData[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(pointInventory_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

## Parameters (to be estimated). /!\ Given on a non standardised scale for dbh
potentialGrowth = log(4.47); # i.e., in average a tree grows 2.47 cm
dbh_slope = 0.11/135;

pr_slope = -0.34;
pr_slope2 = 0.0078;

params = c(potentialGrowth = potentialGrowth, dbh_slope = dbh_slope, pr_slope = pr_slope, pr_slope2 = pr_slope2)

processError = 5.64 # This is a variance: var /!\
measurementError = sqrt(3) # This is a standard deviation: sd /!\ i.e., the variance is 1

## Simulate dbh
set.seed(123)
if (n_indiv != length(children_index))
	stop("Length mismatch between n_indiv and children_index")

latent_dbh = data.table(dbh0 = treeData[parents_index, dbh], dbh1 = numeric(n_indiv), dbh2 = numeric(n_indiv),
	dbh3 = numeric(n_indiv), dbh4 = numeric(n_indiv), dbh5 = numeric(n_indiv))
ls_plot_id = unique(treeData[, pointInventory_id])

clim_ind = indices[type == "parent", index_clim_start]

for (t in 1:5)
{
	precip = (climate[clim_ind, pr] - mu_pr)/sd_pr
	current_col = paste0("dbh", t - 1)
	current_dbh = latent_dbh[, ..current_col]
	setnames(current_dbh, new = "lala")
	current_dbh = current_dbh[, lala]
	latent_dbh[, (paste0("dbh", t)) := current_dbh + rgamma(n_indiv, shape = (growth_fct(current_dbh, precip, params))^2/processError,
		rate = (growth_fct(current_dbh, precip, params))/processError)]
	clim_ind = clim_ind + 1
}

## Replace the observations by the simulated dummy data; do not forget the measurement error!
treeData[parents_index, dbh := rnorm(.N, mean = latent_dbh[, dbh0], sd = measurementError)]
treeData[children_index, dbh := rnorm(.N, mean = latent_dbh[, dbh5], sd = measurementError)]

(checkData = any(treeData[children_index, dbh] - treeData[parents_index, dbh] < 0))
if (checkData)
{
	problems = which(treeData[children_index, dbh] - treeData[parents_index, dbh] < 0)
	print(treeData[c(parents_index[problems[1]], children_index[problems[1]])])
	print(latent_dbh[problems[1]])
}

lm_pa = lm(treeData[parents_index, dbh] ~ latent_dbh[, dbh0])
lm_ch = lm(treeData[children_index, dbh] ~ latent_dbh[, dbh5])

summary(lm_pa) # Should be a slope of 1, an intercept of 0, and a residual of sqrt(3)
summary(lm_ch) # Should be a slope of 1, an intercept of 0, and a residual of sqrt(3)

(sd_dbh = treeData[, sd(dbh)]) # Do not forget to update it to the dummy data!

#### Stan model
## Define stan variables
# Common variables
maxIter = 2000
n_chains = 3

# Initial value for states only
initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]$latent_dbh_parents))

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_precip = climate[, .N], # Dimension of the climate vector
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

	# Explanatory variables
	precip = climate[, pr], # Annual precipitations (sum over 12 months)
	pr_mu = climate_mu_sd[variable == "pr", mu],
	pr_sd = climate_mu_sd[variable == "pr", sd],

	totalTrunkArea = treeData[, totalTrunkArea]
)

## Compile model
model = cmdstan_model("./growth_3rdOption.stan")

## Run model
results = readRDS("./Tilia_platyphyllos/test_procErrorFixed_3rdOption.rds")
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = round(1*maxIter/3), iter_sampling = round(2*maxIter/3), save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 13, adapt_delta = 0.9)

results$save_output_files(dir = paste0("./", species, "/"), basename = "test_procErrorFixed_3rdOption",
	timestamp = FALSE, random = FALSE)
results$save_object(file = paste0("./", species, "/test_procErrorFixed_3rdOption.rds"))

#### Check-up
## Flags
results$cmdstan_diagnose()
results$print()

aa = results$draws("measureError")

# pdf("measureError.pdf")
lazyTrace(sd_dbh*aa, val1 = measurementError, val2 = mean(sd_dbh*aa), label1 = "real", label2 = "est")
# dev.off()
providedArgs = list(val1 = measurementError, val2 = mean(sd_dbh*aa), label1 = "real", label2 = "est")
## Check first five dbh
row = 1:50

estimated_parents = results$draws(paste0("latent_dbh_parents[", row,"]"))
estimated_parents = apply(estimated_parents, 3, mean)
real = unlist(latent_dbh[row, dbh0])
measured = treeData[parents_index[row], dbh]
x = row

ymin = min(sd_dbh*estimated_parents, real, measured)
ymax = max(sd_dbh*estimated_parents, real, measured)

pdf("50shadesOfDBH.pdf")
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(x, sd_dbh*estimated_parents, pch = 19, col = "#34568B", cex = 4,
	xlab = "Index", ylab = "Diameter at breast height (scaled)", ylim = c(ymin, ymax))
points(x, y = measured, pch = 19, col = "#FA7A35", cex = 2)
points(x, real, pch = 19, col = "#CD212A", cex = 1.5)
legend(x = "topleft", legend = c("Estimated", "Measured", "Real"), fill = c("#34568B", "#FA7A35", "#CD212A"),
	box.lwd = 0)
dev.off()

#### Posterior predictive checking: Can the model give rise to new observations that properly resemble the original data?
## Script to generate new data. Note that 'model' is an unnecessary block here as restart from results. gq stands for generated quantities
# The new observations are simulated as follow:
#	1. Draw the vector of parameters theta (which includes the latent states!)
#	2. Generate the parent observation from the corresponding latent dbh according to the model
#	3. Generate the children observation according to the model (i.e., after n years of yearly latent growth)

gq_script = write_stan_file(
"
functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2)
	{
		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_precip; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_precip - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// Explanatory variables
	vector<lower = 0>[n_precip] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector<lower = 0>[n_obs] totalTrunkArea;
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_precip] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) growth
}

generated quantities { // Checked on paper by hand on 4th March 2022
	array [n_obs] real newObservations;
	{ // Variables declared in nested blocks are local variables, not generated quantities, and thus won't be printed.
		real processError = 5.68/sd(Yobs)^2; // processError is a variance
		real current_latent_dbh;
		int growth_counter = 1;
		real m; // mean gamma

		for (i in 1:n_indiv) // Loop over all the individuals
		{
			newObservations[parents_index[i]] = normal_rng(latent_dbh_parents[i], measureError);
			current_latent_dbh = latent_dbh_parents[i];
			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter];
				growth_counter += 1;
			}
			// The last current_latent_dbh corresponds to the child latent dbh
			newObservations[children_index[i]] = normal_rng(current_latent_dbh, measureError);
		}
	}
}
"
)

gq_model = cmdstan_model(gq_script)

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

dim(generate_quantities$draws()) # iter_sampling * n_chains * n_obs

iter_sampling = results$metadata()$iter_sampling

hist(sd_dbh*generate_quantities$draws("newObservations[1]"))
abline(v = treeData[1, dbh], lwd = 4, col = "blue")

