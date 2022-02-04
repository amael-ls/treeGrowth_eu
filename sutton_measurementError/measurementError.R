
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
library(stringi)

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
# Initiate latend_dbh with reasonable value (by default, stan would generate them between 0 and 2---constraint latent_dbh > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh1", "dbh2")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh1 and dbh2")

	dbh1 = providedArgs[["dbh1"]]
	dbh2 = providedArgs[["dbh2"]]
	
	if (length(dbh1) != length(dbh2))
		stop("Arguments should be of the same length")

	return(list(latent_dbh = abs(rnorm(length(dbh1), mean = (dbh1 + dbh2)/2, sd = 5))))
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

## Keep only reliable growth. (This is because the French data I am studying have been corrected)
treeData = treeData[(growth > 0) & (growth <= 10)]
print(paste0("In total, ", treeData[, .N], " individuals are available for the analysis"))

## Estimate the average yearly growth, assuming growth is gamma-distributed
# Common variables
maxIter = 3000
n_chains = 3

# Data
stanData = list(
	# Number of data
	n_trees = treeData[, .N], # Number of trees

	# Data
	growth = treeData[, growth]
)

# Compile model
model_G = cmdstan_model("./estimateGrowth.stan")

# Run model
results = model_G$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE,
	max_treedepth = 14, adapt_delta = 0.8)

# Check results
results$cmdstan_diagnose()

# Get a pointwise value
mu_G = mean(results$draws("mu"))
sigma_G = mean(results$draws("sigma"))

# Plot results
hist(treeData[, growth], probability = TRUE, ylim = c(0, 0.3))
curve(dgamma(x, mu_G^2/sigma_G, mu_G/sigma_G), add = TRUE, col = "#CD212A", lwd = 2)

## Compute the expected dbh if trees grow at average growth
treeData[, expected_dbh_2016 := dbh0_in_mm + years*mu_G]

print(paste("the average growth is", round(treeData[, mean(growth)], digits = 2), "mm/year"))
print(paste("the sd in growth is", round(treeData[, sd(growth)], digits = 2), "mm/year"))
print(paste("the range Δt is", round(treeData[, min(years)], digits = 2), "-", round(treeData[, max(years)], digits = 2)))

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

## Compile model
model = cmdstan_model("./measurementError.stan")

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

print(paste("The average of the estimated measurement error is", round(mean(error), 2), "mm"))

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
