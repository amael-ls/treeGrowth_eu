
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)


#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Create the cluster #! USE A BASH SCRIPT TO RUN THIS R SCRIPT, EXPORT THE ARRAY_ID
array_id = 46 #! REMOVE THIS LINE WHEN DONE WITH TEST. FAGUS_SYLVATICA
array_id = 71 #! REMOVE THIS LINE WHEN DONE WITH TEST. PICEA_ABIES
array_id = 85 #! REMOVE THIS LINE WHEN DONE WITH TEST. PINUS_SYLVESTRIS
array_id = 150 #! REMOVE THIS LINE WHEN DONE WITH TEST. TILIA PLATYPHYLLOS
# array_id = as.integer(Sys.getenv("ARRAY_ID")) # Each array takes care of one species
# print(paste0("array id = ", array_id))

#### Tool function
## Initiate Y_gen with reasonable values (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "n_latentGrowth", "average_yearlyGrowth", "nbYearsGrowth", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh_parents, n_latentGrowth, average_yearlyGrowth, nbYearsGrowth, and normalise")

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	n_latentGrowth = providedArgs[["n_latentGrowth"]]
	average_yearlyGrowth = providedArgs[["average_yearlyGrowth"]]
	nbYearsGrowth = providedArgs[["nbYearsGrowth"]]
	normalise = providedArgs[["normalise"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	if (!normalise)
		Warning("This function has been coded for normalise only. The parameters potential growth, etc are scaled with sd_dbh = 135")

	if (any(average_yearlyGrowth == 0))
	{
		warning("Some average yearly growth were 0. They have been replaced by 0.5")
		average_yearlyGrowth[average_yearlyGrowth == 0] = 0.5
	}

	n_indiv = length(dbh_parents)
	Y_gen = rgamma(n_indiv, dbh_parents^2/0.5, dbh_parents/0.5) # Average = dbh_parents, variance = 0.5

	latent_growth_gen = numeric(n_latentGrowth)
	counter_growth = 0
	for (i in 1:n_indiv) # Not that this forbid trees to shrink
	{
		for (j in 1:nbYearsGrowth[i])
		{
			counter_growth = counter_growth + 1
			latent_growth_gen[counter_growth] = rgamma(n = 1, shape = average_yearlyGrowth[i]^2/0.5, rate = average_yearlyGrowth[i]^2/0.5)
		}
	}

	if (any(latent_growth_gen == 0))
	{
		warning("Some generated latent growth were 0. There have been replaced by 0.1 (before standardising)")
		latent_growth_gen[latent_growth_gen == 0] = 0.1
	}
	
	# All the following initialisations are based on a previous run on Tilia platyphyllos. The set sd is 4 times the estimated sd
	potentialGrowth = rnorm(n = 1, mean = -3.98, sd = 0.08)
	dbh_slope = rnorm(n = 1, mean = 0.10, sd = 0.04)
	
	pr_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)
	pr_slope2 = rnorm(n = 1, mean = 0, sd = 0.04)
	tas_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)
	tas_slope2 = rnorm(n = 1, mean = 0, sd = 0.04)
	competition_slope = rnorm(n = 1, mean = -0.09, sd = 0.04)

	measureError = rgamma(n = 1, shape = 0.011^2/0.0015, rate = 0.011/0.0015)
	processError = rgamma(n = 1, shape = 0.0009^2/0.00012, rate = 0.0009/0.00012)

	# Normalise dbh
	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
		latent_growth_gen = latent_growth_gen/sd_dbh
	}

	return(list(latent_dbh_parents = Y_gen, latent_growth = latent_growth_gen,
		potentialGrowth = potentialGrowth, dbh_slope = dbh_slope,
		pr_slope = pr_slope, pr_slope2 = pr_slope2,
		tas_slope = tas_slope, tas_slope2 = tas_slope2,
		competition_slope = competition_slope))
}

## Function to compute the mean and sd of given variables in a data table
normalisation = function(dt, colnames = names(df), folder = "./", filename = "normalisation.rds", rm_na = TRUE, ...)
{
	if (rm_na)
		print("Warning: rm_na is activated, normalisation won't take NA into account")

	if (!is.data.table(dt))
		stop("This function is written for data table only")

	if (stri_sub(str = folder, from = stri_length(folder)) != "/")
		folder = paste0(folder, "/")

	if (!all(colnames %in% names(dt)))
	{
		warning(paste0("The following columns do not exist in the provided data and are ignored:\n- ",
			paste0(colnames[!(colnames %in% names(dt))], collapse = "\n- ")))
		colnames = colnames[colnames %in% names(dt)]
	}

	providedArgs = list(...)
	providedArgs_names = names(providedArgs)
	if ("indices" %in% providedArgs_names)
	{
		if (!any(c("col_ind", "col_ind_start", "col_ind_end") %in% providedArgs_names))
			stop("Indices provided without any column selected")
		
		ind = providedArgs[["indices"]]
		
		if ("col_ind" %in% providedArgs_names)
		{
			col_ind = providedArgs[["col_ind"]]
			if (!(col_ind %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind))

			rowsToKeep = indices[, ..col_ind]

			if (any(c("col_ind_start", "col_ind_end") %in% providedArgs_names))
				warning("col_ind_start or col_ind_end ignored")
		}

		if (("col_ind_start" %in% providedArgs_names) & !("col_ind" %in% providedArgs_names))
		{
			if (!("col_ind_end" %in% providedArgs_names))
				stop("A starting index is provided but there is no stopping index")
			
			col_ind_start = providedArgs[["col_ind_start"]]
			if (!(col_ind_start %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind_start))

			col_ind_start = providedArgs[["col_ind_start"]]

			col_start = unique(indices[[col_ind_start]])

			col_ind_end = providedArgs[["col_ind_end"]]
			if (!(col_ind_end %in% names(indices)))
				stop(paste("Indices does not contain a column named", col_ind_end))
			
			col_ind_end = providedArgs[["col_ind_end"]]

			col_end = unique(indices[[col_ind_end]])

			if (length(col_start) != length(col_end))
				stop("Starting and ending indices length mismatches")

			rowsToKeep = integer(length = sum(col_end - col_start + 1))
			count = 1
			for (i in 1:length(col_start))
			{
				rowsToKeep[count:(count + col_end[i] - col_start[i])] = col_start[i]:col_end[i]
				count = count + col_end[i] - col_start[i] + 1
			}
		}
	}
	
	n = length(colnames)
	mu_sd = data.table(variable = character(n), mu = numeric(n), sd = numeric(n))

	if (!("indices" %in% providedArgs_names))
		mu_sd[, c("variable", "mu", "sd") := .(colnames, as.matrix(dt[, lapply(.SD, mean, na.rm = rm_na), .SDcols = colnames])[1,],
			as.matrix(dt[, lapply(.SD, sd, na.rm = rm_na), .SDcols = colnames])[1,])]

	if ("indices" %in% providedArgs_names)
	{

		mu_sd[, c("variable", "mu", "sd") := .(colnames, as.matrix(dt[rowsToKeep, lapply(.SD, mean, na.rm = rm_na), .SDcols = colnames])[1,],
			as.matrix(dt[rowsToKeep, lapply(.SD, sd, na.rm = rm_na), .SDcols = colnames])[1,])]
	}

	saveRDS(mu_sd, file = paste0(folder, filename))
	print(paste0("files containing coefficients saved at: ", folder, filename))
}

#### Load data
## Paths
# mainFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
mainFolder = "/home/amael/project_ssm/inventories/FR IFN/processed data/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

# clim_folder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/climateData/Chelsa/yearlyAverage/"
clim_folder = "/home/amael/project_ssm/climateData/Chelsa/yearlyAverage/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
# treeData[, .N, by = speciesName_sci][1:20,]
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((array_id < 1) | (array_id > length(ls_species)))
	stop(paste0("Array id = ", array_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

species = ls_species[array_id]
print(paste("Script running for species:", species))
treeData = treeData[speciesName_sci == species]
savingPath = paste0("./", species, "/")
if (!dir.exists(savingPath))
	dir.create(savingPath)

## Climate
climate = readRDS(paste0(clim_folder, "FR_reshaped_climate.rds"))

## Set-up indices
indices = readRDS(paste0(mainFolder, species, "_indices.rds"))

# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1], type := "parent"]

# Compute the number of growing years per individual
indices[, nbYearsGrowth := max(year) - min(year), by = .(pointInventory_id, tree_id)]

checkUp = all(indices[, nbYearsGrowth == index_clim_end - index_clim_start])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

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
print(paste("Number of latent dbh:", n_hiddenState))
n_latentGrowth = n_hiddenState - n_indiv

parents_index = treeData[, .I[which.min(year)], by = .(pointInventory_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(pointInventory_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(pointInventory_id, tree_id)][, V1]
# not_parent_index = 1:indices[.N, index_gen]
# not_parent_index = not_parent_index[!(not_parent_index %in% indices[type == "parent", index_gen])]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimension mismatch between children_index and number of children")

# if (length(not_parent_index) != indices[.N, index_gen] - n_indiv)
# 	stop("Dimension mismatch between not_parent_index, n_hiddenState, and n_indiv")

#### Compute and save normalising constantes
normalisation(dt = treeData, colnames = "dbh", folder = savingPath, filename = "dbh_normalisation.rds")
normalisation(dt = climate, colnames = c("pr", "tas"), folder = savingPath, filename = "climate_normalisation.rds", indices = indices,
	col_ind_start = "index_clim_start", col_ind_end = "index_clim_end")

climate_mu_sd = readRDS(paste0(savingPath, "climate_normalisation.rds"))

#### Stan model
## Define stan variables
# Common variables
maxIter = 2100
n_chains = 3

# Initial values for states only
average_yearlyGrowth = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_yearlyGrowth) != n_indiv)
	stop("Dimensions mismatch between average_yearlyGrowth and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
	n_latentGrowth = n_latentGrowth, average_yearlyGrowth = average_yearlyGrowth, nbYearsGrowth = nbYearsGrowth,
	normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]][["latent_growth"]]))

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_climate = climate[, .N], # Dimension of the climate vector
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

	tas = climate[, tas], # Annual average temperature (average over 12 months)
	tas_mu = climate_mu_sd[variable == "tas", mu],
	tas_sd = climate_mu_sd[variable == "tas", sd],

	totalTreeWeight = treeData[parents_index, totalTreeWeight] # Includes also the weight of other species!
)

saveRDS(object = stanData, file = paste0(savingPath, "stanData.rds"))

## Compile model
model = cmdstan_model("./growth_3rdOption.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = round(2*maxIter/3), iter_sampling = round(maxIter/3), save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 15, adapt_delta = 0.9)

time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-", time_ended, ".rds"))

results$cmdstan_diagnose()

results$print(c("potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"competition_slope", "measureError", "processError"))



#* --------------------------------------------------------------------------------------------
#! ********************************************************************************************
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================    IT WOULD BE GOOD THAT THE SCRIPT STOPS HERE    ================ ODOT#
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! ********************************************************************************************
#* --------------------------------------------------------------------------------------------


lazyTrace(results$draws("measureError"))
lazyTrace(results$draws("processError"))

#### Compute residuals: compare data versus latent states with obs error
n_rep = results$metadata()$iter_sampling  * results$num_chains()

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

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dt_dharma[seq(1, .N, by = n_rep), rep_dbh/sd(treeData$dbh)]) # treeData[, dbh/sd(dbh)]

# pdf("child_residuals.pdf", height = 6, width = 9)
plot(forDharma)
dev.off()


par(mfrow = c(2,1))
lazyTrace(latent_dbh_array[, , "latent_dbh[1]"], val1 = treeData[1, dbh]/sd(treeData[, dbh]), main = "tree 1, t0", scaling = sd(treeData$dbh))
lazyTrace(latent_dbh_array[, , "latent_dbh[6]"], val1 = treeData[2, dbh]/sd(treeData[, dbh]), main = "tree 1, t1", scaling = sd(treeData$dbh))
dev.off()

par(mfrow = c(2,1))
lazyTrace(latent_dbh_array[, , "latent_dbh[7]"], val1 = treeData[3, dbh]/sd(treeData[, dbh]), main = "tree 2, t0", scaling = sd(treeData$dbh))
lazyTrace(latent_dbh_array[, , "latent_dbh[12]"], val1 = treeData[4, dbh]/sd(treeData[, dbh]), main = "tree 2, t1", scaling = sd(treeData$dbh))
dev.off()

####! CRASH TEST ZONE

## Function to compute growth, the data table must be sorted by year within tree id and plot id
computeGrowth = function(dt, col = "growth", byCols = c("pointInventory_id", "tree_id"))
{
	while (col %in% names(dt))
	{
		newcol = paste0(col, rnorm(1))
		warning(paste0("The name `", col, "` is already used in the data table. The result was stored in col `", newcol, "` instead"))
		col = newcol
	}
	dt[, (col) := (shift(dbh, n = 1, type = "lead", fill = NA) - dbh)/(shift(year, n = 1, type = "lead", fill = NA) - year), by = byCols]

	if (!("deltaYear" %in% names(dt)))
		dt[, deltaYear := shift(year, n = 1, type = "lead", fill = NA) - year, by = byCols]
}

computeGrowth(treeData)
dataGrowth = na.omit(treeData)[, .(year, tree_id, pointInventory_id, dbh, growth, deltaYear)]

pr_avg = climate[, mean(pr), by = pointInventory_id]
setnames(pr_avg, "V1", "pr")

dataGrowth = pr_avg[dataGrowth, on = "pointInventory_id"]

sd_dbh = sd(treeData[, dbh])
mu_pr = climate_mu_sd[variable == "pr", mu]
sd_pr = climate_mu_sd[variable == "pr", sd]

dataGrowth[, pr_norm := (pr - mu_pr)/sd_pr]

params = results$draws(c("potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2"))
params = apply(params, 3, mean)

growth_fct = function(dbh0, params, precip, scaling_dbh)
{
	potentialGrowth = params["potentialGrowth"]
	dbh_slope = params["dbh_slope"]
	pr_slope = params["pr_slope"]
	pr_slope2 = params["pr_slope2"]

	return(scaling_dbh*exp(potentialGrowth + dbh_slope*dbh0/scaling_dbh + pr_slope*precip + pr_slope2*precip^2))
}

plot(dataGrowth[, dbh], dataGrowth[, growth], pch = 19)
curve(growth_fct(x, params, 0, sd_dbh), col = "#178F92", lwd = 3, from = 0, to = 700, add = TRUE)

lazyTrace(results$draws("potentialGrowth"))
lazyTrace(results$draws("dbh_slope"))

latent_1_6 = results$draws(paste0("latent_dbh[", 1:6, "]"))
latent_1_6 = apply(latent_1_6, 3, mean)
x = 2000:2005

ymin = min(sd_dbh*latent_1_6, treeData[1, dbh], treeData[2, dbh])
ymax = max(sd_dbh*latent_1_6, treeData[1, dbh], treeData[2, dbh])

lm1 = lm(sd_dbh*latent_1_6 ~ x)
lm2 = lm(c(treeData[1, dbh], treeData[2, dbh]) ~ c(min(x), max(x)))

op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(2000:2005, sd_dbh*latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)", ylim = c(ymin, ymax))
points(x = 2000, y = treeData[1, dbh], pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = treeData[2, dbh], pch = 19, col = "#CD212A", cex = 2)
abline(lm1, col = "#EC9706", lwd = 2)
abline(lm2, col = "#178F92", lwd = 2)
dev.off()


lazyTrace(results$draws("alpha"))
lazyTrace(results$draws("beta"))

var_params = results$draws(c("alpha","beta"))
var_params = apply(var_params, 3, mean)

var_fct = function(dbh, params, scaling)
{
	alpha = params["alpha"]
	beta = params["beta"]
	return (scaling^2*exp(alpha + beta * dbh/scaling))
}

curve(var_fct(x, var_params, sd_dbh), from = 0, to = 600)

####! END CRASH TEST ZONE


#? -------------------------------------------------------------------------------------------
#* ######################    Restart run from an already tuned model    ######################
#? -------------------------------------------------------------------------------------------

#### FROM https://discourse.mc-stan.org/t/saving-reusing-adaptation-in-cmdstanr/19166/44
get_inits = function(chain_id){
	warmup_draws = results$draws(inc_warmup = TRUE)
	final_warmup_value = warmup_draws[maxIter/2, chain_id, 2:(dim(warmup_draws)[3])]
	(
		final_warmup_value
		%>% tibble::as_tibble(.name_repair = function(names){
				dimnames(final_warmup_value)$variable
			})
		%>% tidyr::pivot_longer(cols=dplyr::everything())
		%>% tidyr::separate(name, into = c('variable','index'), sep = "\\[", fill = 'right')
		%>% dplyr::group_split(variable)
		%>% purrr::map(
			.f = function(x){
				(
					x
					%>% dplyr::mutate(
						index = stringr::str_replace(index,']','')
					)
					%>% tidyr::separate(
						index
						, into = c('first','second')
						, sep = ","
						, fill = 'right'
						, convert = T
					)
					%>% dplyr::arrange(first,second)
					%>% (function(x){
						out = list()
						if(all(is.na(x$second))){
							out[[1]] = x$value
						}else{
							out[[1]] = matrix(
								x$value
								, nrow = max(x$first)
								, ncol = max(x$second)
							)
						}
						names(out) = x$variable[1]
						return(out)
					})
				)
			}
		)
		%>% unlist(recursive = FALSE)
		%>% return()
	)
}

results2 = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = 0, iter_sampling = maxIter,
	adapt_engaged = FALSE,
	inv_metric = results$inv_metric(matrix = FALSE),
	step_size = results$metadata()$step_size_adaptation,
	init = get_inits)

results2$cmdstan_diagnose()

latent_1_6 = getParams(results, paste0("latent_dbh[", 1:6, "]"))
x = 2000:2005
min_y = min(135.137*latent_1_6, treeData[1:2, dbh])
max_y = max(135.137*latent_1_6, treeData[1:2, dbh])

pdf("./latent_real_2ndOption.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(2000:2005, 135.137*latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)", ylim = c(min_y, max_y))
points(x = 2000, y = treeData[1, dbh], pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = treeData[2, dbh], pch = 19, col = "#CD212A", cex = 2)
dev.off()

lazyTrace(results$draws("latent_dbh[1]"))

plot_title = ggplot2::ggtitle("latent_dbh[1]")
bayesplot::mcmc_trace(135.137*results$draws("latent_dbh[1]")) + plot_title

plot_title = ggplot2::ggtitle("Traces for measureError")
bayesplot::mcmc_trace(results$draws("measureError")) + plot_title


####! CRASH TEST ZONE 2
climate[, avg_pr := mean(pr), by = c("pointInventory_id")]
climate[, avg_tas := mean(tas), by = c("pointInventory_id")]
climate[, avg_tasmin := mean(tasmin), by = c("pointInventory_id")]
climate[, avg_tasmax := mean(tasmax), by = c("pointInventory_id")]

aa = climate[treeData, on = c("pointInventory_id", "year")]
aa[seq(2, .N, by = 2), G := dbh - aa[seq(1, .N, by = 2), dbh]]
aa[seq(2, .N, by = 2), firstTotTrunkA := aa[seq(1, .N, by = 2), totalTrunkArea]]
aa = na.omit(aa)

aa[, avg_pr_Z := (avg_pr - mean(avg_pr))/sd(avg_pr)]
aa[, avg_tas_Z := (avg_tas - mean(avg_tas))/sd(avg_tas)]
aa[, avg_tasmin_Z := (avg_tasmin - mean(avg_tasmin))/sd(avg_tasmin)]
aa[, avg_tasmax_Z := (avg_tasmax - mean(avg_tasmax))/sd(avg_tasmax)]

pr = lm(aa[, G] ~ aa[, avg_pr_Z] + aa[, avg_pr_Z^2])
tas = lm(aa[, G] ~ aa[, avg_tas])
tasmin = lm(aa[, G] ~ aa[, avg_tasmin])
tasmax = lm(aa[, G] ~ aa[, avg_tasmax])
comp = lm(aa[, G] ~ aa[, firstTotTrunkA])
comp_log10 = lm(aa[, log10(G + 10)] ~ aa[, firstTotTrunkA])

together = lm(aa[, G] ~ aa[, avg_pr_Z] + aa[, avg_tas_Z])

plot(aa[, firstTotTrunkA], aa[, log10(G + 10)], pch = 19, cex = 0.5, xlab = "Trunk area 5 years ago", ylab = "5 years dbh increment")
abline(comp_log10, col = "orange", lwd = 2)

par(mfrow = c(2,2))
plot(comp_log10)
dev.off()

pdf("avg_growth-versus-avg_precip.pdf", width = 6, height = 6)
plot(aa[, avg_pr], aa[, G], pch = 19, cex = 0.5, xlab = "Average precip over 5 years", ylab = "5 years dbh increment")
abline(pr, col = "orange", lwd = 2)
dev.off()

plot(aa[, avg_tas], aa[, G], pch = 19, cex = 0.5, xlab = "Average precip over 5 years", ylab = "5 years dbh increment")
abline(tas, col = "orange", lwd = 2)

plot(aa[, avg_tasmin], aa[, G], pch = 19, cex = 0.5, xlab = "Average precip over 5 years", ylab = "5 years dbh increment")
abline(tasmin, col = "orange", lwd = 2)

plot(aa[, avg_tasmax], aa[, G], pch = 19, cex = 0.5, xlab = "Average precip over 5 years", ylab = "5 years dbh increment")
abline(tasmax, col = "orange", lwd = 2)

library(mgcv)
pr_gam = gam(G ~ s(avg_pr_Z, bs="cr"), data = aa)

# Note that mod_gam2$model is the data that was used in the modeling process, 
# so it will have NAs removed.
testdata = data.frame(avg_pr_Z = seq(-2, 3, length = 300))
fits = predict(pr_gam, newdata = testdata, type='response', se = TRUE)
predicts = data.frame(testdata, fits) %>% 
	dplyr::mutate(lower = fit - 1.96*se.fit, upper = fit + 1.96*se.fit)

plot_mod_gam2_response = ggplot(aes(x=avg_pr_Z,y=fit), data=predicts) +
	geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
	geom_line(color='#00aaff', lwd = 2) +
	geom_point(aes(x=avg_pr_Z,y=G), data = aa, cex = 0.5) +
	theme_classic()

# TEMPERATURE
tas_gam = gam(G ~ s(avg_tas_Z, bs="cr"), data = aa)

# Note that mod_gam2$model is the data that was used in the modeling process, 
# so it will have NAs removed.
testdata = data.frame(avg_tas_Z = seq(-4, 4, length = 400))
fits = predict(tas_gam, newdata = testdata, type='response', se = TRUE)
predicts = data.frame(testdata, fits) %>% 
	dplyr::mutate(lower = fit - 1.96*se.fit, upper = fit + 1.96*se.fit)

plot_mod_gam2_response = ggplot(aes(x=avg_tas_Z,y=fit), data=predicts) +
	geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
	geom_line(color='#00aaff', lwd = 2) +
	geom_point(aes(x=avg_tas_Z,y=G), data = aa, cex = 0.5) +
	theme_classic()

####! END CRASH TEST ZONE 2
