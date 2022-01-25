
#### Aim of prog: Fit the French growth data. (This objective will be updated to part of Europe later)
#! IMPORTANT REMARK:
# There are two ways to use Y_generated_0:
#	1. It is a data, and therefore can be used at any time to help convergence on either the parents only or on the hidden states too
#		depending on the size of the vector and the likelihood
#	2. It is an initial condition, which help to start the hidden states. This also means that Y_generated_0 "is forgotten" as iterations
#		goes by
#
# cd projects/def-dgravel/amael/postdoc/bayForDemo/growth/
#
#! NOTES
# The file growth-2022-01-18_21h19.rds I HAVE NO IDEA WHAT WORKS. This model did not converge, but the algo is correctly tuned. So try to run it longer! What are the priors?
# The file with measureError = 3.0/135.137 is growth-2022-01-19_20h26.rds. This model converged.
# The file with measureError = 1.0/135.137 is growth-2022-01-23_21h01.rds. This model converged.

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(doParallel)
library(cmdstanr)
library(magrittr)
library(stringi)

if (!("callr" %in% installed.packages()))
	install.packages("callr", repos = "http://cran.us.r-project.org")

if (!("future" %in% installed.packages()))
	install.packages("future", repos = "http://cran.us.r-project.org")

#### Create the cluster
## Cluster variables
array_id = 150 #! REMOVE THIS LINE WHEN DONE WITH TEST
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) # Each array takes care of one year and process all the variables for that year
print(paste0("array id = ", array_id))

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
			Y_gen[count + j] = Y_gen[count + j - 1] + average_G[i] + rgamma(1, shape = 2.5, rate = 5) # mean = 0.5, var = 0.1

		count = count + years_indiv[i];
	}

	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
	}
	return(list(latent_dbh = Y_gen))
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
mainFolder = "Tilia_platyphyllos/" # Folder to use when running on local computer
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

# clim_folder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/climateData/Chelsa/yearlyAverage/"
clim_folder = mainFolder # Folder to use when running on local computer
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((array_id < 1) | (array_id > length(ls_species)))
	stop(paste0("Array id = ", array_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

species = ls_species[array_id]
treeData = treeData[speciesName_sci == species]
savingPath = paste0("./", species, "/")
if (!dir.exists(savingPath))
	dir.create(savingPath)

## Climate
climate = readRDS(paste0(clim_folder, "FR_reshaped_climate.rds"))

## inidices
indices = readRDS(paste0(mainFolder, species, "_indices.rds"))

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

#### Compute and save normalising constantes
normalisation(dt = treeData, colnames = "dbh", folder = savingPath, filename = "dbh_normalisation.rds")
normalisation(dt = climate, colnames = "pr", folder = savingPath, filename = "climate_normalisation.rds", indices = indices,
	col_ind_start = "index_clim_start", col_ind_end = "index_clim_end")

climate_mu_sd = readRDS(paste0(savingPath, "climate_normalisation.rds"))

#### Stan model
## Define stan variables
# Common variables
maxIter = 1.5e3
n_chains = 3

# Initial value for states only
average_G = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen], normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]))

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

# jpeg("./initVal.jpg", height = 1080, width = 1080, quality = 100)
# plot(1:indices[.N, index_gen], initVal_Y_gen[[1]]$latentState, pch = 19, col = "#34568B",
# 	xlab = "Tree index", ylab = "Diameter at breast height (in mm)")
# points(x = indices[type == "parent", index_gen], y = treeData[parents_index, dbh], pch = 19, col = "#FA7A35")
# points(x = indices[type == "child", index_gen], y = treeData[children_index, dbh], pch = 19, col = "#CD212A")
# dev.off()

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

	# tas = climate[, tas], # Annual average temperature (average over 12 months)
	# tas_mu = climate_mu_sd[variable == "tas", mu],
	# tas_sd = climate_mu_sd[variable == "tas", sd],

	totalTrunkArea = treeData[, totalTrunkArea],

	# Diffuse initialisation for the parents
	initialParents = initVal_Y_gen[[1]]$latent_dbh[indices[type == "parent", index_gen]]

	# Parameter for parallel calculus
	# grainsize = 1
)

## Compile model
# model = cmdstan_model("growth.stan")
model = cmdstan_model("./growth_2ndOption.stan")

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 12, adapt_delta = 0.9)

time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-", time_ended, ".rds"))

diagnose = results$cmdstan_diagnose()

# stanfit <- rstan::read_stan_csv(fit$output_files())
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


####! CRASH
climate[, avg_pr := mean(pr), by = c("pointInventory_id")]
climate[, avg_tas := mean(tas), by = c("pointInventory_id")]
climate[, avg_tasmin := mean(tasmin), by = c("pointInventory_id")]
climate[, avg_tasmax := mean(tasmax), by = c("pointInventory_id")]

aa = climate[treeData, on = c("pointInventory_id", "year")]
aa[seq(2, .N, by = 2), G := dbh - aa[seq(1, .N, by = 2), dbh]]
aa = na.omit(aa)

aa[, avg_pr_Z := (avg_pr - mean(avg_pr))/sd(avg_pr)]
aa[, avg_tas_Z := (avg_tas - mean(avg_tas))/sd(avg_tas)]
aa[, avg_tasmin_Z := (avg_tasmin - mean(avg_tasmin))/sd(avg_tasmin)]
aa[, avg_tasmax_Z := (avg_tasmax - mean(avg_tasmax))/sd(avg_tasmax)]

pr = lm(aa[, G] ~ aa[, avg_pr_Z] + aa[, avg_pr_Z^2])
tas = lm(aa[, G] ~ aa[, avg_tas])
tasmin = lm(aa[, G] ~ aa[, avg_tasmin])
tasmax = lm(aa[, G] ~ aa[, avg_tasmax])



allTogether = lm(aa[, G] ~ aa[, avg_pr_Z] + aa[, avg_tas_Z])

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
