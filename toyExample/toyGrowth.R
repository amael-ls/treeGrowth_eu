# library(data.table)
# # Dummy data
# dt = data.table(year = c(2000, 2001, 2003, 2001, 2005, 2000, 2008),
# 	g1 = c(1, 1, 1, 2, 2, 3, 3), g2 = c(88, 88, 88, 88, 88, 54, 54))

# # Set up new col to foo
# dt[, newCol := "foo"]

# # Correct the value for the minimal year, by group g1 and g2
# dt[dt[, .I[which.min(year)], by = .(g1, g2)][, V1], newCol := "bar"]

# dt[, newCol := c("foo", "bar")[1 + (year == min(year))], .(g1, g2)]

# dt[year == min(year), newCol := 'bar' ,.(g1, g2)] # Does not work

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(doParallel)
library(rstan)

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

#### Create the cluster
## Cluster variables
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

print(paste("number of cores:", future::availableCores()))

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
		Y_gen[count + 1] = rgamma(1, shape = dbh_parents[i]^2, rate = dbh_parents[i]) # Mean = dbh_parents[i], Variance = 1
		for (j in 2:years_indiv[i])
			Y_gen[count + j] = Y_gen[count + j - 1] + average_G[i] + rgamma(1, shape = 0.25, rate = 0.5)

		count = count + years_indiv[i];
	}

	return(list(Y_generated = Y_gen))
}

#### Load data
# treeData = readRDS("./toyData.rds")
# precip = readRDS("./precipitation.rds")
# indices = readRDS("./indices.rds")

treeData = readRDS("./toyData2.rds")
precip = readRDS("./precipitation2.rds")
indices = readRDS("./indices2.rds")

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

n_indiv = unique(treeData[, .(tree_id, pointInventory_id)])[, .N]
print(paste("Number of individuals:", treeData[, .N]))

nbYearsPerIndiv = unique(indices[, .(tree_id, plot_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
if (length(nbYearsPerIndiv) != n_indiv)
	stop("Dimension mismatch between nbYearsPerIndiv and n_indiv")

parents_index = treeData[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(pointInventory_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(pointInventory_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

#### Stan model
## Define stan variables
# Common variables
maxIter = 2e3
n_chains = 4

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_precip = length(precip), # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_hiddenState = indices[.N, index_gen], # Dimension of the state space vector
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsPerIndiv = nbYearsPerIndiv, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in Y_generated
	climate_index = indices[type == "parent", index_precip_start], # Index of the climate associated to each parent

	# Observations
	Yobs = treeData[, dbh],

	# Explanatory variable
	precip = precip # Precipitations
)

# Initial value for states only
average_G = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen])

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

## Compile stan model
model = stan_model(file = "./toy.stan")

## Run model
start = proc.time()

results = stan(file = "toy.stan", data = stanData, cores = n_chains,
	iter = maxIter, chains = n_chains, init = initVal_Y_gen)

proc.time() - start

saveRDS(results, "test_2.rds")

stopCluster(cl)


#! CRASH TEST ZONE
results = readRDS("./test_2.rds")

results_summary = summary(results)
saveRDS(results_summary, "test_2_summary.rds")

results_summary_merged = results_summary$summary
head(results_summary_merged)
class(results_summary_merged)

dt = as.data.table(results_summary_merged, keep.rownames = TRUE)
dt_1.1 = dt[Rhat > 1.1]
dt_1.8 = dt[Rhat > 1.8]

dt_1.1[, .N]
dt_1.8[, .N]

pdf("test_2-intercept.pdf")
traceplot(results, pars = "intercepts[1]", inc_warmup = TRUE)
dev.off()

pdf("test_2-slopes_precip.pdf")
traceplot(results, pars = "slopes_precip[1]", inc_warmup = TRUE)
dev.off()

pdf("test_2-slopes_dbh.pdf")
traceplot(results, pars = "slopes_dbh[1]")
dev.off()

pdf("test_2-Y_generated.pdf")
traceplot(results, pars = paste0("Y_generated[", 1:4,"]"))
dev.off()

pdf("test_2-intercepts_sd.pdf")
traceplot(results, pars = "intercepts_sd")
dev.off()

pdf("test_2-precip_sd.pdf")
traceplot(results, pars = "precip_sd")
dev.off()

pdf("test_2-quad_slopes_precip[14897].pdf")
traceplot(results, pars = "quad_slopes_precip[14897]")
dev.off()


# A = init_fun(chain_id = 1, dbh_parents = treeData[parents_index, dbh],
# 	years_indiv = nbYearsPerIndiv, average_G = average_G,
# 	n_hiddenState = indices[.N, index_gen])

# A = A[[1]]
# A[1:6]
# treeData[1:2, .(year, dbh)]

# A[7:12]
# treeData[3:4, .(year, dbh)]

# graphics.off()
# plot(A[1:5], A[2:6], xlim = range(A), ylim = range(A), pch = 20)
# for (i in 2:n_indiv)
# {
# 	points(A[(6*(i - 1) + 1):(6*i - 1)], A[(6*(i - 1) + 2):(6*i)], xlim = range(A), ylim = range(A), pch = 20)
# 	if (i %% 1000 == 0)
# 		print(paste0(i*100/n_indiv, " % done"))
# }

#! END CRASH TEST ZONE