
#### Aim of prog: Investigate on the measurement error
## Explanations
# Around 950 trees have been measured by two people in Sutton Canada. It was done this way:
#	1. The first person measure many trees (more than a hundred in a row) while the second person takes notes
#	2. Then, the roles are reversed: person 2 measures and person 1 takes notes.
# There should have enough measurements in a row to make sure that the second person is not influenced by the first
# 
## Names column data:
#	- Arbre = tree id number
#	- Esp = species
#	- Etat = state (V = alive, there are only living trees)
#	- Multi = is it a multi-trunk or not. N for No
#	- The last two columns are the values collected by the two measurers. Note that the names of the measurers appear in the data provided
#		by the lab Gravel (Sherbrooke, Canada). To preserve the anonymity, these names were replaced by dbh1_in_mm, dbh2_in_mm
#	- The two persons measured the diameters of trees in mm at the height 1.37m (breast height).

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(bayesplot)
library(cmdstanr)
library(stringi)

#### Tool functions
## Initiate latend_dbh with reasonable value (by default, stan would generate them between 0 and 2---constraint latent_dbh > 0)
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

#### Load data
treeData = readRDS("./trees_remeasured.rds")

treeData[, mean(dbh1_in_mm)]
treeData[, sd(dbh1_in_mm)]

treeData[, mean(dbh2_in_mm)]
treeData[, sd(dbh2_in_mm)]

plot(treeData[, dbh1_in_mm], treeData[, dbh2_in_mm])

#### Estimate parameters
## Prepare stan data
# Common variables
maxIter = 3000
n_chains = 3

# Initialialisation
latent_dbh_gen = lapply(1:n_chains, init_fun, dbh1 = treeData[, dbh1_in_mm], dbh2 = treeData[, dbh2_in_mm])
maxLatent_dbh = ceiling(max(unlist(latent_dbh_gen)))

length(latent_dbh_gen)

for (i in 1:n_chains)
	print(range(latent_dbh_gen[[i]]))

# Data top provide
stanData = list(
	# Number of data
	n_trees = treeData[, .N], # Number of trees

	# Data
	dbh1 = treeData[, dbh1_in_mm],
	dbh2 = treeData[, dbh2_in_mm]
	# maxDBH_allowed = max(max(treeData[, dbh1_in_mm]), max(treeData[, dbh2_in_mm]), maxLatent_dbh) + 5 # The +5 is to make it 'loose'
)

## Compile model
model = cmdstan_model("./measurementError.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE,
	max_treedepth = 10, adapt_delta = 0.8)

## Check-up
results$cmdstan_diagnose()

## Saving
time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = "./", basename = paste0("sutton-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0("./", "sutton-", time_ended, ".rds"))

#### Get latent_dbh and error
mcmc_trace(results$draws("latent_dbh[1]"))
lazyTrace(results$draws("latent_dbh[2]"), val1 = treeData[2, dbh1_in_mm], val2 = treeData[2, dbh2_in_mm])
results$print()
results$print(variables = "latent_dbh")

error = results$draws("error")

####! CRASH TEST ZONE
dbetabinom = function(x, size, shape1, shape2, log = FALSE)
{
	lpmf = lchoose(size, x) +
		lbeta(x + shape1, size - x + shape2) -
		lbeta(shape1, shape2)
	if (log)
		return (lpmf);
		
	return (exp(lpmf));
}

rbetabinom = function(n, size, shape1, shape2)
	rbinom(n, size = size, prob = rbeta(n, shape1 = shape1, shape2 = shape2))

y = dbetabinom(x = 1:10, size = 5, shape1 = 600, shape2 = 400)
plot(1:10, y, type = "l", lwd = 2, col = "blue")
points(1:10, y, pch = 19)

meanCompute = function(n, shape1, shape2)
	return (n * shape1/(shape1 + shape2))

varCompute = function(n, shape1, shape2)
	return (n*shape1*shape2*(shape1 + shape2 + n)/((shape1 + shape2)^2*(shape1 + shape2 + 1)));

shape1Compute = function(n, mean, var)
	return ((-mean^3 + mean^2*n - mean*var)/(mean^2 - mean*n + var*n));

shape2Compute = function(n, mean, var)
	return ((mean - n)*(mean^2 - mean*n + var)/(mean^2 - mean*n + var*n));

aa = rbetabinom(1e3, 10, 600, 400)

mean(aa)
var(aa)
meanCompute(10, 600, 400)
varCompute(10, 600, 400)

shape1Compute(10, mean(aa), var(aa)) # Should be 600
shape1Compute(10, 6, 2.421578) # Should be 600
shape1Compute(10, 5.966, 2.349193) # -225.3685, really far from 600! This is based on a sample of 1e3
shape1Compute(10, 5.98846, 2.418951) # 776.7303, really far from 600, but positive! This is based on a sample of 1e5!
shape2Compute(10, mean(aa), var(aa)) # Should be 400

meanCompute(10, shape1Compute(10, mean(aa), var(aa)), shape2Compute(10, mean(aa), var(aa)))
varCompute(10, shape1Compute(10, mean(aa), var(aa)), shape2Compute(10, mean(aa), var(aa)))

curve(meanCompute(x, 600, 400), 0, 300)
curve(meanCompute(x, 30000, 0.8), 0, 300, col = "blue", add = TRUE)
curve(varCompute(x, 30000, 0.8), 0, 30)
curve(varCompute(x, 30000, 0.8), 0, 30, col = "red", add = TRUE)

# The two parameters shape1 and shape2 are highly sensitive to mean and var. This prevent me to use the beta binomial distribution

####! END CRASH TEST ZONE
