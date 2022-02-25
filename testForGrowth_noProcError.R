
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
		stop("You must provide chain_id, Yobs, years_indiv, and normalise")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	normalise = providedArgs[["normalise"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	n_indiv = length(dbh_parents)

	Y_gen = rgamma(n_indiv, shape = dbh_parents^2, rate = dbh_parents) # Mean = dbh_parents, Variance = 1

	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
	}
	return(list(latent_dbh = Y_gen))
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
	potentialGrowth = params[["potentialGrowth"]]
	dbh_slope = params[["dbh_slope"]]
	pr_slope = params[["pr_slope"]]
	pr_slope2 = params[["pr_slope2"]]
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

## Parameters (that to be estimated). /!\ Given on a non standardised scale for dbh
potentialGrowth = log(2.47); # i.e., in average a tree grows 2.47 cm
dbh_slope = 0.11/135;

pr_slope = -0.34;
pr_slope2 = 0.0078;

params = c(potentialGrowth = potentialGrowth, dbh_slope = dbh_slope, pr_slope = pr_slope, pr_slope2 = pr_slope2)

processError = 5.68 # This is a variance: var /!\
measurementError = sqrt(3) # This is a standard deviation: sd /!\ i.e., the variance is 3

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
	latent_dbh[, (paste0("dbh", t)) := current_dbh + growth_fct(current_dbh, precip, params)]
	clim_ind = clim_ind + 1
}

## Replace the observations by the simulated dummy data; do not forget the measurement error!
treeData[parents_index, dbh := rnorm(.N, mean = latent_dbh[, dbh0], sd = measurementError)]
treeData[children_index, dbh := rnorm(.N, mean = latent_dbh[, dbh5], sd = measurementError)]

sd_dbh = treeData[, sd(dbh)]

#### Stan model
## Define stan variables
# Common variables
maxIter = 1800
n_chains = 3

# Initial value for states only
initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
	normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

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
model = cmdstan_model("./growth_noProcError.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 13, adapt_delta = 0.9)

results$save_output_files(dir = paste0("./", species, "/"), basename = "test_noProcError",
	timestamp = FALSE, random = FALSE)
results$save_object(file = paste0("./", species, "/test_noProcError.rds"))

#### Check-up
## Flags
results$cmdstan_diagnose()
results$print()

## Check first five dbh
row = 1:5

estimated_parents = results$draws(paste0("latent_dbh_parents[", row,"]"))
estimated_parents = apply(estimated_parents, 3, mean)
real = unlist(latent_dbh[row, dbh0])
measured = treeData[parents_index[row], dbh]
x = row

ymin = min(sd_dbh*estimated_parents, real, measured)
ymax = max(sd_dbh*estimated_parents, real, measured)

op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(x, sd_dbh*estimated_parents, pch = 19, col = "#34568B", cex = 4,
	xlab = "Index", ylab = "Diameter at breast height (scaled)", ylim = c(ymin, ymax))
points(x, y = measured, pch = 19, col = "#FA7A35", cex = 2)
points(x, real, pch = 19, col = "#CD212A", cex = 1.5)
legend(x = "topleft", legend = c("Estimated", "Measured", "Real"), fill = c("#34568B", "#FA7A35", "#CD212A"),
	box.lwd = 0)

#### Compute residuals: compare data versus latent states with obs error
n_rep = results$metadata()$iter_sampling * results$num_chains()

dt_dharma = data.table(
	rep_latent_id = rep(1:n_indiv, each = n_rep),
	rep_dbh = rep(treeData[parents_index, dbh], each = n_rep),
	sampled = numeric(n_rep * treeData[, .N]/2))

latent_dbh_array = results$draws("latent_dbh_parents") # dimension: iter_sampling * n_chains * number latent states (= n_indiv)
dt_dharma[, sampled := myPredictObs(latent_dbh_array, rep_latent_id, regex = "latent_dbh_parents"), by = rep_latent_id]

sims = matrix(data = dt_dharma[, sampled], nrow = n_rep, ncol = treeData[, .N]/2) # each column is for one data point
sims = t(sims) # Transpose the matrix for dharma
dim(sims)

if (ncol(sims) != n_rep)
	stop("Dimensions mismatch")

dt_dharma[, diff := rep_dbh/sd_dbh - sampled]
meanDiff = dt_dharma[, mean(diff), by = rep_latent_id]
print(paste0("The mean is ", round(meanDiff[, mean(V1)], digits = 4)))
plot(meanDiff[, rep_latent_id], meanDiff[, V1], xlab = "Latent tree id", ylab = "Average difference: measured - simulated", pch = 19)

hist(meanDiff[, V1])
hist(dt_dharma[, diff])

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dt_dharma[seq(1, .N, by = n_rep), rep_dbh/sd(treeData$dbh)]) # treeData[, dbh/sd(dbh)]

pdf("residuals_parent_test2.pdf")
plot(forDharma)
dev.off()

#### Retrieve the parameters back
# If there is no mistake (see notebook), then the relation between some stan parameters and the parameters for non standardised dbh is:
#	potentialGrowth_stan = potentialGrowth_R - log(sd(dbh))
#	slope_dbh_stan = sd(dbh)*slope_dbh_R,
# Where *_stan are for the stan estimates, while *_R are for the real values in R
# Note that there is no transformation required for the climate parameters, only the dbh_slope and the intercept are influenced

(potentialGrowth_stan = mean(results$draws("potentialGrowth")))
potentialGrowth - log(sd(treeData[, dbh]))

(dbh_slope_stan = mean(results$draws("dbh_slope")))
sd(treeData[, dbh])*dbh_slope

(pr_slope_stan = mean(results$draws("pr_slope")))
pr_slope

(pr_slope2_stan = mean(results$draws("pr_slope2")))
pr_slope2

lazyTrace(results$draws("potentialGrowth"), val1 = potentialGrowth - log(sd(treeData[, dbh])))
lazyTrace(results$draws("dbh_slope"), val1 = sd(treeData[, dbh])*dbh_slope)
lazyTrace(results$draws("pr_slope"), val1 = pr_slope)
lazyTrace(results$draws("pr_slope2"), val1 = pr_slope2)

## Log-likelihood function for test 1
loglik = function(data, latent_dbh_parents, parents_index, children_index, nbYearsPerIndiv, params, normalised_precip, climate_index)
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
	
	ll = sum(dnorm(data[parents_index, dbh], mean = latent_dbh_parents, sd = sqrt(3), log = TRUE)) +
		sum(dnorm(data[children_index, dbh], mean = computed_latent_child, sd = sqrt(3), log = TRUE)) +
		dnorm(params["potentialGrowth"], mean = 0, sd = 100, log = TRUE) + dnorm(params["dbh_slope"], mean = 0, sd = 5, log = TRUE) +
		sum(dnorm(pr_slope, mean = 0, sd = 5, log = TRUE)) + sum(dnorm(pr_slope2, mean = 0, sd = 5, log = TRUE))
	return (ll)
}

estimated_parents = results$draws("latent_dbh_parents")
estimated_parents = apply(estimated_parents, 3, mean)

vec_pG = seq(potentialGrowth - 0.1, potentialGrowth + 0.1, length.out = 40)
vec_dbhSlope = seq(dbh_slope - 0.01/135, dbh_slope + 0.01/135, length.out = 40)

ll = matrix(data = 0, nrow = length(vec_pG), ncol = length(vec_dbhSlope))

dim(ll)
r = 1
for (pG in vec_pG)
{
	c = 1
	for (dbhSlope in vec_dbhSlope)
	{
		currentParams = c(potentialGrowth = pG, dbh_slope = dbhSlope, pr_slope = pr_slope, pr_slope2 = pr_slope2)
		ll[r, c] = loglik(data = treeData, latent_dbh_parents = sd_dbh*estimated_parents, parents_index = parents_index,
		children_index = children_index, nbYearsPerIndiv = nbYearsPerIndiv,
		params = c(potentialGrowth = pG, dbh_slope = dbhSlope, pr_slope = pr_slope_stan, pr_slope2 = pr_slope2_stan),
		normalised_precip = (climate[, pr] - mu_pr)/sd_pr, climate_index = indices[type == "parent", index_clim_start])
		c = c + 1
	}
	if (r %% 5 == 0)
		print(r)
	r = r + 1
}

x_sol = dbh_slope_stan/sd_dbh
y_sol = potentialGrowth_stan + log(sd_dbh)
x_real = params["dbh_slope"]
y_real = params["potentialGrowth"]

pdf(paste0("./", species, "/test_noProcError.pdf"))
filled.contour(x = vec_dbhSlope, y = vec_pG, z = t(ll), xlab = "dbh slope", ylab = "potential growth",
	plot.axes = {
		axis(1); axis(2);
		points(x_sol, y_sol, col = "#1D1F54", pch = 19, cex = 2); text(x_sol, y_sol, pos = 2, labels = "Estimated")
		points(x_real, y_real, col = "#1D7FDF", pch = 19, cex = 2); text(x_real, y_real, pos = 4, labels = "Real solution")
	})
dev.off()

#### Check if the stan program can generate correctly the children latent states using the block generated_quantities
# More information can be found at https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html

## Define the script to generate the latent states. Note that model is actually an unnecessary block here. gq stands for generated quantities
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
		int<lower = 2> n_hiddenState; // Dimension of the state space vector
		int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
		int<lower = 2, upper = n_obs> nbYearsPerIndiv[n_indiv]; // Number of years for each individual

		// Indices
		int<lower = 1, upper = n_obs - 1> parents_index[n_indiv]; // Index of each parent in the observed data
		int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
		int<lower = 1, upper = n_precip - 1> climate_index[n_indiv]; // Index of the climate associated to each parent

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

		// real<lower = 0> processError; // Constrained by default, realistically not too small
		// real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

		vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	}

	generated quantities {
		real computed_latent[n_hiddenState];
		{ // Variables declared in nested blocks are local variables, not generated quantities, and thus won't be printed.
			int count = 0;
			
			for (i in 1:n_indiv) // Loop over all the individuals
			{
				count = count + 1;
				computed_latent[count] = latent_dbh_parents[i];
				for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
				{
					count = count + 1;
					computed_latent[count] = computed_latent[count - 1] + growth(computed_latent[count - 1],
						normalised_precip[climate_index[i] + j - 2], potentialGrowth, dbh_slope, pr_slope, pr_slope2);
				}
			}
		}
	}
	"
)

gq_model = cmdstan_model(gq_script)

generate_quantities = gq_model$generate_quantities(results, data = stanData, parallel_chains = n_chains)

dim(generate_quantities$draws()) # maxIter * n_chains * n_hiddenStates

sampling_start = results$metadata()$iter_warmup + 1
n_iter = sampling_start + results$metadata()$iter_sampling - 1
gq_dbh5 = generate_quantities$draws()[sampling_start:n_iter, , indices[type == "child", index_gen]]
dim(gq_dbh5) # iter_sampling * n_chains * n_indiv

n_rep = results$metadata()$iter_sampling * results$num_chains()

dt = data.table(rep_id_latent = rep(indices[type == "child", index_gen], each = n_rep),
	observed = rep(treeData[children_index, dbh], each = n_rep),
	latent = rep(latent_dbh[, dbh5], each = n_rep),
	sampled = numeric(n_indiv*n_rep))

dt[, sampled := myPredictObs(draws_array = gq_dbh5, id_latent = rep_id_latent, regex = "computed_latent"), by = rep_id_latent]
dt[, sampled := sampled*sd_dbh]

dt[, res_obs := sampled - observed]
dt[, res_lat := latent - observed]


pdf("residuals_hist_test2.pdf")
hist(dt[, res_obs])
dev.off()

sd(dt[, res_obs])
sqrt(3)

sims = matrix(data = dt[, sampled], nrow = n_rep, ncol = treeData[, .N]/2) # each column is for one data point
sims = t(sims) # Transpose the matrix for dharma
dim(sims)

if (ncol(sims) != n_rep)
	stop("Dimensions mismatch")

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dt[seq(1, .N, by = n_rep), observed]) # treeData[, dbh/sd(dbh)]

pdf("residuals_child_test2.pdf")
plot(forDharma)
dev.off()

