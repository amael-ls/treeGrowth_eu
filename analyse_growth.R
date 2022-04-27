
#### Aim of prog: Analysing results (check-up residuals, plots)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(reticulate)
library(cmdstanr)
library(stringi)
library(png)

use_python("/usr/bin/python3")

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
species = "Tilia platyphyllos" # "Tilia platyphyllos", "Fagus sylvatica"
path = paste0("./", species, "/")
run = 1
info_lastRun = getLastRun(path = path, run = run)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]

isDBH_normalised = TRUE
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

## Compile simulation generator
gq_model = cmdstan_model("./generate_posteriorSimulations.stan")

## Generate simulations
# Access data
n_chains = results$num_chains()
iter_sampling = results$metadata()$iter_sampling
n_obs = stanData$n_obs
n_hiddenState = stanData$n_latentGrowth + stanData$n_indiv
indices = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "indices.rds"))

stanData$nfi_id = unique(indices[, .(plot_id, tree_id, nfi_index)])[, nfi_index]

if (length(stanData$nfi_id) != stanData$n_indiv)
	stop("Dimensions mismatch")

# for (i in 1:length(stanData))
# 	if (sum(is.na(stanData[[i]])) != 0)
# 		print(i) #! THERE ARE NA IN THE pH!!!!!

# names(stanData)
# stanData = stanData[-27]

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
dim(newObservations_array) # iter_sampling * n_chains * n_obs
sum(is.na(newObservations_array))

latent_dbh_array = generate_quantities$draws("latent_dbh_parentsChildren")
dim(latent_dbh_array) # iter_sampling * n_chains * n_obs
sum(is.na(latent_dbh_array))

dt_dharma[, simulated_observations := reshapeDraws(newObservations_array, rep_id, regex = "newObservations"), by = rep_id]
dt_dharma[, simulated_observations := sd_dbh*simulated_observations]

dt_dharma[, latent_dbh := reshapeDraws(latent_dbh_array, rep_id, regex = "latent_dbh_parentsChildren"), by = rep_id]
dt_dharma[, latent_dbh := sd_dbh*latent_dbh]

dt_dharma[, residuals_obs := rep_dbh - simulated_observations]

dt_dharma[rep_id %in% stanData$parents_index, type := "parents"]
dt_dharma[rep_id %in% stanData$children_index, type := "children"]

setorderv(x = dt_dharma, cols = c("type", "rep_id"), order = c(-1, 1)) # The -1 is to have parents first

dt_dharma[, mean(residuals_obs), by = type]
print(paste("The residuals' variance is", round(var(dt_dharma[, residuals_obs]), 3)))

# jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_scatter.jpg"), quality = 50)
# plot(dt_dharma[, rep_id], dt_dharma[, residuals_obs], pch = '.', col = "#A1A1A122")
# abline(v = stanData$n_indiv, lwd = 2, col = "#CD212A")
# axis(3, at = n_obs/4, "Parents", las = 1)
# axis(3, at = 3*n_obs/4, "Children", las = 1)
# dev.off()

mpl = import("matplotlib")
mpl$use("Agg") # Stable non interactive back-end
plt = import("matplotlib.pyplot")
mpl$rcParams['agg.path.chunksize'] = 0 # Disable error check on too many points

plt$figure()
plt$plot(dt_dharma[, rep_id], dt_dharma[, residuals_obs], '.', c = "#A1A1A122", markersize = 1)
plt$axvline(x = stanData$n_indiv, linewidth = 2, color = "#CD212A")
plt$savefig(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_scatter.png"))
plt$close(plt$gcf())

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_hist.pdf"))
hist(dt_dharma[, residuals_obs])
dev.off()

qq = qqnorm(dt_dharma[, residuals_obs], pch = 1, frame = FALSE, plot.it = FALSE)

jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "qqplot.jpg"))
plot(qq)
qqline(dt_dharma[, residuals_obs], col = "#34568B", lwd = 3)
dev.off()

## Check that the process error does not show any pattern with any predictor
latentG_residuals_array = generate_quantities$draws("latentG_residuals") #! Used to be procError_sim
mean(latentG_residuals_array)

dim(latentG_residuals_array) # iter_sampling * n_chains * n_latentGrowth

latentG_residuals_vec = as.vector(latentG_residuals_array) # The order is array[,1,1], array[,2,1], ..., array[,n_chain,1], array[,2,1], ...
length(latentG_residuals_vec)

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "latentG_residuals_check_hist.pdf"))
hist(latentG_residuals_vec)
dev.off()

latentG_residuals_avg = apply(latentG_residuals_array, 3, mean) # Average latentG_residuals for each latent growth (or dbh)
length(latentG_residuals_avg)

index_notLastMeasure = 1:n_hiddenState
index_notLastMeasure = index_notLastMeasure[!(index_notLastMeasure %in% indices[stanData$last_child_index, index_gen])]

if (length(index_notLastMeasure) != length(latentG_residuals_avg))
	stop("Dimension mismatch between latentG_residuals_avg and index_notLastMeasure")

latent_dbh = apply(generate_quantities$draws("yearly_latent_dbh"), 3, mean)
latent_dbh = latent_dbh[index_notLastMeasure]

jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "latentG_residuals_vs_dbh_check.jpg"), quality = 50)
plot(latent_dbh, latentG_residuals_avg, type = "l")
dev.off()

cor(latent_dbh, latentG_residuals_avg)

# #### Crash test zone
# dt_dharma[which.min(dt_dharma[, residuals_obs])]
# mainFolder = "/home/amael/project_ssm/inventories/growth/"
# treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))
# treeData = treeData[speciesName_sci == species]
# treeData[(75.7577 < dbh) & (dbh < 75.7578), .(plot_id, tree_id, year)]

# treeData[plot_id %in% c("france_1327882", "france_743874") & tree_id %in% c(9, 4)]

# 3: Tilia platyphyllos FR IFN  france_743874       4  2007 85.92677 2.492116 45.32078      21.073138  france wfo-0000456948
# 4: Tilia platyphyllos FR IFN  france_743874       4  2012 89.12677 2.492116 45.32078      22.557006  france wfo-0000456948
# 5: Tilia platyphyllos FR IFN  france_743874       4  2017 75.75775 2.492116 45.32078      22.931993  france wfo-0000456948
