
#### Aim of prog: Investigate on the measurement error
## Explanations
# The data used in this script are from the lab Integrative Ecology (PI: Dominique GRAVEL, https://ielab.recherche.usherbrooke.ca)
#
# 945 trees have been measured by two people in Sutton Canada. It was done this way:
#	1. The first person measure many trees (more than a hundred in a row) while the second person takes notes
#	2. Then, the roles are reversed: person 2 measures and person 1 takes notes.
# There should have enough measurements in a row to make sure that the second person is not influenced by the first
# This was done the last week of May 2016
# I also got data from previous campaign in Sutton from 2011 to 2013, and from 2011 to 2016
# 
## Names column data:
#	- Arbre = tree id number
#	- Esp = species
#	- Etat = state (V = alive, there are only living trees)
#	- Multi = is it a multi-trunk or not. N for No
#	- The last two columns are the values collected by the two measurers. Note that the names of the measurers appeared in the data provided
#		by the lab Gravel (Sherbrooke, Canada). To preserve the anonymity, these names were replaced by dbh1_in_mm, dbh2_in_mm
#	- The two persons measured the diameters of trees in mm at the height 1.37m (breast height).
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(bayesplot)
library(cmdstanr)
library(MetBrewer)
library(stringi)
library(DHARMa)

#### Tool functions
## Related to checking/cleaning data
# To check that each tree has only one species associated (i.e., individuals do not change species)
checkSpecies = function(dt, usedKey = "Arbre", option = "simplify")
{
	if ("usedKey" %in% names(dt))
		warning("Data table might be confused by the argument usedKey. It is both an argument of the function and a column name")

	dt[, species_nb := length(unique(Esp)), by = usedKey]
	if (dt[species_nb > 1, .N] != 0)
	{
		warning("Some trees have more than one species!")
		if (option == "simplify")
		{
			warning("You chose the option to simplify, species will be rewritten as unknown")
			dt[species_nb > 1, Esp := "unknown"]
		}
	}
	dt[, species_nb := NULL]
}

# To check that each tree id is unique (not that multi-trunk trees can have more than one row. Use the ... argument with toKeep)
checkIndividuals = function(dt, usedKey = "Arbre", option = "remove", ...)
{
	if ("usedKey" %in% names(dt))
		warning("Data table might be confused by the argument usedKey. It is both an argument of the function and a column name")
	
	dt[, nb_indiv := .N, by = usedKey]
	if (dt[, max(nb_indiv) != 1])
	{
		warning("Some individuals have more than one row of data. First simplification with function `unique'")
		dt = unique(dt)
		dt[, nb_indiv := .N, by = usedKey]
		
		if (dt[, max(nb_indiv) != 1])
		{
			warning("Despite applying unique, there are still some individuals with at least two rows of data")
			if (option == "remove")
			{
				warning("You chose the option to remove the concerned rows (except those that are to kept, set-up by the key)")
				dt = dt[nb_indiv == 1]
			}
		}
	}
	dt[, nb_indiv := NULL]
	return (dt) # No other choice than returning because I am removing some data. Not possible to do it by reference
}

## Related to Stan
## Function to create the individual indices for the state-space model (stan)
indices_state_space = function(trees_NFI)
{
	start = 0
	end = 0

	trees_NFI[, dbh0_ind := 0]
	trees_NFI[, dbh2016_ind := 0]

	for (i in 1:trees_NFI[, .N])
	{
		start = end + 1
		end = start + trees_NFI[i, round_years]
		trees_NFI[i, c("dbh0_ind", "dbh2016_ind") := .(start, end)]
	}
	if (trees_NFI[.N, dbh2016_ind] != sum(trees_NFI[, round_years + 1]))
		stop("Indices were miscomputed")
}

# Initiate latend_dbh with reasonable value (by default, stan would generate them between 0 and 2---constraint latent_dbh > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh0", "average_G", "sd_G", "n_indiv", "n_states", "parents", "children")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh0 and average_G, sd_G, n_states, parents, and children")

	dbh0 = providedArgs[["dbh0"]]
	average_G = providedArgs[["average_G"]]
	sd_G = providedArgs[["sd_G"]]
	n_indiv = providedArgs[["n_indiv"]]
	n_states = providedArgs[["n_states"]]
	parents = providedArgs[["parents"]]
	children = providedArgs[["children"]]

	if (length(dbh0) != n_indiv)
		stop("Error in length dbh0")

	if (length(parents) != n_indiv)
		stop("Error in length parents")

	if (length(children) != n_indiv)
		stop("Error in length children")

	latent_dbh = numeric(length = n_states)

	latent_dbh[parents] = dbh0
	for (indiv in 1:n_indiv)
		for (ind in (parents[indiv] + 1):(children[indiv]))
			latent_dbh[ind] = latent_dbh[ind - 1] + rgamma(1, shape = average_G^2/sd_G, rate = average_G/sd_G)

	return(list(latent_dbh = latent_dbh))
}

# To plot traces of mcmcs from cmdstan fit objects
lazyTrace = function(draws, ...)
{
	if (!all(class(draws) %in% c("draws_array", "draws", "array")))
		stop("The class of draws should be compatible with stan (draws_array, draws, array)")
	
	n_chains = dim(draws)[2]
	n_iter = dim(draws)[1]
	cols = c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")

	min_val = min(draws)
	max_val = max(draws)
	plot(0, pch = "", xlim = c(0, n_iter), ylim = c(min_val, max_val), axes = TRUE,
		xlab = "", ylab = "", bg = "transparent")

	for (chain in 1:n_chains)
		lines(1:n_iter, draws[, chain,], type = "l", col = cols[chain])	
	
	providedArgs = list(...)
	nbArgs = length(providedArgs)

	if (nbArgs != 0)
	{
		ls_names = names(providedArgs)
		val_ind = stri_detect(str = ls_names, regex = "val")

		if (any(val_ind))
		{
			for (val in ls_names[val_ind])
				abline(h = providedArgs[[val]], col = "#CD212A")
		}
	}
}

#### Loading and cleaning data
## Load three dataset corresponding to trees measured twice the same day in May 2016, and measures before and after 2016
treeData_remeasured = readRDS("./trees_remeasured.rds")
print(paste("The data set contains", treeData_remeasured[, .N], "measures"))
print(paste("The data set contains", length(treeData_remeasured[, unique(Arbre)]), "individuals"))

treeData_campaign_1 = readRDS("./trees_campaign1.rds")
setnames(treeData_campaign_1, "DHP_mm_", "dbh_in_mm")

treeData_campaign_2 = readRDS("./trees_campaign2.rds")
setnames(treeData_campaign_2, "DHP1.mm.", "dbh_in_mm")

## Checking the data and correcting them
# Error in the dates the 19th May 2016 for the dataset treeData_campaign_2
treeData_campaign_2[Date == "206-05-19", Date := "2016-05-19"] # I know from the lab Gravel that 2016 is the only possible date
treeData_campaign_2[Date == "20165-19", Date := "2016-05-19"]

# Checking species
checkSpecies(treeData_remeasured)
checkSpecies(treeData_campaign_1)
checkSpecies(treeData_campaign_2)

# Checking uniqueness of the individuals (see crash test zone to keep few more data). I will remove the multi-trunk trees
multiTrunk = treeData_remeasured[Multi == "O", Arbre] # List the trees that have multi-trunk! Those trees for sure have more than one line

treeData_remeasured = treeData_remeasured[!(Arbre %in% multiTrunk)]
treeData_campaign_1 = treeData_campaign_1[!(Arbre %in% multiTrunk)]
treeData_campaign_2 = treeData_campaign_2[!(Arbre %in% multiTrunk)]

treeData_remeasured = checkIndividuals(treeData_remeasured, usedKey = "Arbre")
treeData_campaign_1 = checkIndividuals(treeData_campaign_1, usedKey = "Arbre")
treeData_campaign_2 = checkIndividuals(treeData_campaign_2, usedKey = "Arbre")

## Keep only living trees
treeData_remeasured = treeData_remeasured[Etat == "V"]
treeData_campaign_1 = treeData_campaign_1[Etat == "V"]
treeData_campaign_2 = treeData_campaign_2[Etat == "V"]

## Remove rows containing empty fields for treeData_campaign_2
treeData_campaign_2[, dbh_in_mm := as.integer(dbh_in_mm)]
treeData_campaign_2 = na.omit(treeData_campaign_2)

## Last check-up
if (treeData_remeasured[, length(unique(Arbre))] != treeData_remeasured[, .N])
	warning("Check the uniqueness of the data again!")

if (treeData_campaign_1[, length(unique(Arbre))] != treeData_campaign_1[, .N])
	warning("Check the uniqueness of the data again!")

if (treeData_campaign_2[, length(unique(Arbre))] != treeData_campaign_2[, .N])
	warning("Check the uniqueness of the data again!")

## Few descriptions
treeData_campaign_1[, range(Date)]
treeData_campaign_2[, range(Date)]

## Keep only the years for the dates
# Add a date to treeData_remeasured (done the last week of May 2016). I chose the 25th May 2016 (middle of the week)
treeData_remeasured[, Date := as.Date("2016-05-25", format = "%Y-%m-%d")]
treeData_campaign_1[, Date := as.Date(Date, format = "%Y-%m-%d")]
treeData_campaign_2[, Date := as.Date(Date, format = "%Y-%m-%d")]

#### Merging the dataset
## Set key "Arbre" for merging. This key is unique since I removed the multi-trunk trees
setkey(treeData_remeasured, "Arbre")
setkey(treeData_campaign_1, "Arbre")
setkey(treeData_campaign_2, "Arbre")

## Keep only the column of interest
treeData_remeasured = treeData_remeasured[, .(Arbre, Esp, dbh1_in_mm, dbh2_in_mm, Date)]
treeData_campaign_1 = treeData_campaign_1[, .(Arbre, dbh_in_mm, Date)]
treeData_campaign_2 = treeData_campaign_2[, .(Arbre, dbh_in_mm, Date)]

## Merging data sets
# Actually, I just discovered that the second dataset is useless... The interesting trees are the same as in the first dataset!
treeData = treeData_campaign_1[treeData_remeasured, on = "Arbre"] 

## Remove NA
treeData = na.omit(treeData)

treeData[, mean(dbh1_in_mm)]
treeData[, sd(dbh1_in_mm)]

treeData[, mean(dbh2_in_mm)]
treeData[, sd(dbh2_in_mm)]

plot(treeData[, dbh1_in_mm], treeData[, dbh2_in_mm])

#### Compute the growth beween the first measure and the remeasurements
## Change colnames
setnames(treeData, c("Arbre", "dbh_in_mm", "Date", "i.Date"), c("tree_id", "dbh0_in_mm", "date_begin", "date_end"))

## Compute time diff
treeData[, years := as.numeric(date_end - date_begin)/365.25]

## Data yearly growth
treeData[, growth := ((dbh1_in_mm + dbh2_in_mm)/2 - dbh0_in_mm)/years]

## Keep only reliable growth, that is to say positive growth but
treeData_growth = treeData[(growth > 0) & (growth <= 10)]
scaling_dbh_G = sd(treeData_growth[, dbh0_in_mm])
print(paste0("In total, ", treeData_growth[, .N], " individuals are available for the estimation of growth"))

## Estimate the average yearly growth, assuming growth is gamma-distributed
# Common variables
set.seed(1969-08-18) # Woodstock seed
maxIter = 3000
n_chains = 3

# Data
stanData = list(
	# Number of data
	n_trees = treeData_growth[, .N], # Number of trees

	# Data
	growth = treeData_growth[, growth],
	dbh0 = treeData_growth[, dbh0_in_mm]
)

# Compile model
model_G = cmdstan_model("./estimateGrowth.stan")

# Run model
results_G = model_G$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE,
	max_treedepth = 14, adapt_delta = 0.9)

# Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results_G$save_output_files(dir = "./", basename = paste0("model_growth-", time_ended), timestamp = FALSE, random = TRUE)
results_G$save_object(file = paste0("./", "model_growth-", time_ended, ".rds"))

## Check results
# Is there any stan flags
results_G$cmdstan_diagnose()
results_G$print()

# Check residuals, ptw = pointwise estimation, mean in our case
a_ptw = mean(results_G$draws("a"))
b_ptw = mean(results_G$draws("b"))
c_ptw = mean(results_G$draws("c"))

g_ptw = mean(results_G$draws("g"))
h_ptw = mean(results_G$draws("h"))
i_ptw = mean(results_G$draws("i"))

# mu_ptw = mean(results_G$draws("mu"))
# sigma_ptw = mean(results_G$draws("sigma"))
n_rep = 3000
dt_dharma = data.table(rep_dbh = rep(treeData_growth[, dbh0_in_mm/sd(dbh0_in_mm)], each = n_rep),
	sampled = numeric(n_rep * treeData_growth[, .N]))

mu_G = function(x, g, h, i, scaling_sd)
	return(g*exp(-exp(h - i*x/scaling_sd)))

var_G = function(x, a, b, c, scaling_sd)
	return(exp(a*(x/scaling_sd)^2 + b*x/scaling_sd + c))

dt_dharma[, sampled := rgamma(.N,
	shape = mu_G(rep_dbh, g_ptw, h_ptw, i_ptw, 1)^2/var_G(rep_dbh, a_ptw, b_ptw, c_ptw, 1),
	rate = mu_G(rep_dbh, g_ptw, h_ptw, i_ptw, 1)/var_G(rep_dbh, a_ptw, b_ptw, c_ptw, 1))]

sims = matrix(data = dt_dharma[, sampled], nrow = n_rep, ncol = treeData_growth[, .N]) # each column is for one data point
sims = t(sims) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = treeData_growth[, growth])

pdf("suttonGrowth_residuals.pdf", height = 6, width = 9)
plot(forDharma)
dev.off()

# Plot growth in function of diameter
pdf("suttonGrowth.pdf", width = 6, height = 6)
plot(treeData_growth[, dbh0_in_mm], treeData_growth[, growth], xlab = "dbh (in mm)", ylab = "Growth (in mm/yr)",
	pch = 19)
curve(mu_G(x, g = g_ptw, h = h_ptw, i = i_ptw, scaling_sd = treeData_growth[, sd(dbh0_in_mm)]), add = TRUE, col = "#E8731E", lwd = 4)
dev.off()

## Rounding the number of growth years and get few informations on the data
treeData[, round_years := round(years)]

print(paste("the average growth over one year is",
	round(mean(mu_G(treeData[, dbh0_in_mm], g_ptw, h_ptw, i_ptw, treeData[, sd(dbh0_in_mm)])), digits = 2),
	"mm/year"))
print(paste("the sd in growth is", round(treeData[, sd(growth)], digits = 2), "mm/year"))
print(paste("the range Δt is", round(treeData[, min(years)], digits = 2), "-", round(treeData[, max(years)], digits = 2)))

#? #######################################################################################################
#* ##################      ESTIMATION OF THE MEASUREMENT ERROR FOR THE FRENCH DATA      ##################
#? #######################################################################################################
#### Explanation:
# The French data have been corrected (by the agency in charge of measuring the trees). Therefore, the measurement error correspond to how the
# tape is placed, and which tension is put into the tape. The obvious errors (such as typos that create insane growth) have been corrected.
# Therefore, we filter the data from Sutton to remove insane growth and estimate the remaining variability in the measurements.

treeData_french = treeData[(growth > 0) & (growth <= 10)]
indices_state_space(treeData_french)

# hist(treeData_growth[, dbh0_in_mm], breaks = seq(50, 630, by = 10), freq = FALSE)
# ll = mean(treeData_growth[, dbh0_in_mm])
# qq = var(treeData_growth[, dbh0_in_mm])
# uu = log(ll^2/sqrt(qq + ll^2))
# ss = sqrt(log(10*v/m^2 + 1)) 
# curve(dgamma(x, ll^2/qq, ll/qq), add = TRUE, col = "red", lwd = 4)
# curve(dlnorm(x, uu, ss), add = TRUE, lwd = 4, col = "black")

#### Estimate parameters
## Prepare stan data
# Common variables
maxIter = 3000
n_chains = 3

# Initialialisation
latent_dbh_gen = lapply(1:n_chains, init_fun, dbh0 = treeData_french[, dbh0_in_mm], average_G = 2.0, sd_G = 3.0,
	n_states = sum(treeData_french[, round_years + 1]), n_indiv = length(treeData_french[, unique(tree_id)]),
	parents = treeData_french[, dbh0_ind], children = treeData_french[, dbh2016_ind])

length(latent_dbh_gen)

for (i in 1:n_chains)
	print(range(latent_dbh_gen[[i]]))

# Data top provide
stanData = list(
	# Number of data
	n_trees = treeData_french[, .N], # Number of trees
	n_states = sum(treeData_french[, round_years + 1]), # Number of states

	# indices
	index_parents = treeData_french[, dbh0_ind], # index of the measurement corresponding to the first campaign (2011-2013)
	index_children = treeData_french[, dbh2016_ind], # index of the measurement corresponding to the campaign in 2016

	# Parameters growth
	a = a_ptw, # For var_G
	b = b_ptw, # For var_G
	c = c_ptw, # For var_G

	g = g_ptw, # For mu_G
	h = h_ptw, # For mu_G
	i = i_ptw, # For mu_G

	scaling = scaling_dbh_G, # Scaling dbh that was USED to parameterise growth

	# Data
	dbh0 = treeData_french[, dbh0_in_mm], # dbh measured in the first campaign (2011-2013) neither by person 1 nor 2
	dbh1 = treeData_french[, dbh1_in_mm], # dbh measured by person 1 in 2016
	dbh2 = treeData_french[, dbh2_in_mm], # dbh measured by person 2 in 2016
	years = treeData_french[, round_years] # Number of years spent between the first campaign and 2016
)

## Compile model
model = cmdstan_model("./measurementError.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE, init = latent_dbh_gen,
	max_treedepth = 14, adapt_delta = 0.8)

## Check-up
results$cmdstan_diagnose()
results$print()

## Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = "./", basename = paste0("model_error_correctedData-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0("./", "model_error_correctedData-", time_ended, ".rds"))

#### Get latent_dbh and error
lazyTrace(results$draws("latent_dbh[1]"), val1 = treeData_french[1, dbh0_in_mm])
lazyTrace(results$draws("latent_dbh[6]"), val1 = treeData_french[1, dbh1_in_mm], val2 = treeData_french[1, dbh2_in_mm])

error = results$draws("error")
lazyTrace(results$draws("error"))

print(paste("The average of the estimated measurement error (expressed as an sd) is", round(mean(error), 2), "mm"))
print(paste("The average of the estimated measurement error (expressed as a var) is", round(mean(error)^2, 2), "mm"))

####! CRASH TEST ZONE, WHAT FOLLOW IS NOT YET READY
#? ##########################################################################################################
#* ##################      ESTIMATION OF THE MEASUREMENT ERROR FOR NON CORRECTED DATA      ##################
#? ##########################################################################################################
#### Explanation:
# This is for non corrected data, that is to say: I keep 'insane' growth

#### Estimate parameters
## Prepare stan data
# Common variables
maxIter = 3000
n_chains = 3

# Initialialisation
latent_dbh_gen = lapply(1:n_chains, init_fun, dbh1 = treeData[, dbh1_in_mm], dbh2 = treeData[, dbh2_in_mm])

length(latent_dbh_gen)

for (i in 1:n_chains)
	print(range(latent_dbh_gen[[i]]))

# Data top provide
stanData = list(
	# Number of data
	n_trees = treeData[, .N], # Number of trees

	# Data
	dbh1 = treeData[, dbh1_in_mm],
	dbh2 = treeData[, dbh2_in_mm],
	expected_dbh_2016 = treeData[, expected_dbh_2016],
	sd_growth = treeData[, years*sqrt(sigma_G)] # Reminder var[aX] = a^2 var[X], so that sd[aX] = a sd[X]
)

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE,
	max_treedepth = 14, adapt_delta = 0.8)

## Check-up
results$cmdstan_diagnose()

## Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = "./", basename = paste0("sutton-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0("./", "sutton-", time_ended, ".rds"))

#### Get latent_dbh and error
lazyTrace(results$draws("latent_dbh[1]"), val1 = treeData[1, dbh1_in_mm], val2 = treeData[1, dbh2_in_mm], val3 = treeData[1, expected_dbh_2016])
lazyTrace(results$draws("latent_dbh[2]"), val1 = treeData[2, dbh1_in_mm], val2 = treeData[2, dbh2_in_mm], val3 = treeData[2, expected_dbh_2016])
results$print()
results$print(variables = "latent_dbh")

error = results$draws("error")
mcmc_trace(results$draws("error"))
lazyTrace(results$draws("error"))

print(paste("The average of the estimated measurement error (expressed as an sd) is", round(mean(error), 2), "mm"))
print(paste("The average of the estimated measurement error (expressed as a var) is", round(mean(error)^2, 2), "mm"))






####! CRASH TEST ZONE
# The two parameters shape1 and shape2 are highly sensitive to mean and var. This prevent me to use the beta binomial distribution

##? To keep more data (see the cleaining part of the data)

# multiTrunk = treeData_remeasured[Multi == "O", Arbre] # Keep the trees that have multi-trunk! Those trees for sure have more than one line
# treeData_remeasured[, unique_id := as.character(Arbre)]
# treeData_remeasured[Arbre %in% multiTrunk, unique_id := paste0(unique_id, "_", 1:.N), by = Arbre]

# treeData_campaign_1[, unique_id := as.character(Arbre)]
# treeData_campaign_1[Arbre %in% multiTrunk, unique_id := paste0(unique_id, "_", 1:.N), by = Arbre]

# treeData_campaign_2[, unique_id := as.character(Arbre)]
# treeData_campaign_2[Arbre %in% multiTrunk, unique_id := paste0(unique_id, "_", 1:.N), by = Arbre]

# treeData_remeasured = checkIndividuals(treeData_remeasured, usedKey = "unique_id")
# treeData_campaign_1 = checkIndividuals(treeData_campaign_1, usedKey = "unique_id")
# treeData_campaign_2 = checkIndividuals(treeData_campaign_2, usedKey = "unique_id")

# idToKeep = treeData_campaign_1[Arbre %in% multiTrunk]
# idToKeep = treeData_remeasured[idToKeep, on = "Arbre", nomatch = 0]
# idToKeep[, diff := mean(c(dbh1_in_mm, dbh2_in_mm)) - dbh_in_mm, by = i.unique_id]

# union = ?? # Should do it manually since there is not that much trees in idToKeep

#treeData_campaign_1[treeData_campaign_2, on = "Arbre", nomatch = 0]

## Stupid question about rules on probabilities...
aa = runif(n = 1000, min = 0, max = 200)
mu = 2
sigma = 3

dt = data.table(c0 = aa, c1 = numeric(length(aa)), c2 = numeric(length(aa)), c3 = numeric(length(aa)),
	c4 = numeric(length(aa)), c5 = numeric(length(aa)))

for (i in 1:5)
{
	current = paste0("c", i)
	prev = paste0("c", i - 1)

	dt[, (current) := get(prev) + rgamma(.N, shape = mu^2/sigma, rate = mu/sigma)]
}

hist(dt[, c5], probability = TRUE)
bb = rgamma(1e5, (5*mu + aa)^2/(25*sigma), (5*mu + aa)/(25*sigma))
dd = density(bb)
lines(dd)

mean(dt[, c5])
mean(bb)

var(dt[, c5])
var(bb)

# So it seems that my intuition is good. After n iteration of the Markov Chain, the moments are:
#		mean(dbh_n) = n*mu + dbh0
#		sd(dbh_n) = n^2*sigma

####! END CRASH TEST ZONE
