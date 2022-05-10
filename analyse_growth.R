
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

#### Tool functions
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
getLastRun = function(path, begin = "^growth-", extension = ".rds$", format = "ymd", run = NULL, getAll = FALSE, hour = TRUE)
{
	if (format != "ymd")
		stop("Format is not recognised. For now only year-month-day alias ymd works")
	
	if (!is.null(run))
	{
		print(paste("Searching among runs =", run))
		begin = paste0(begin, "run=", run, "-")
	}
	
	ls_files = list.files(path = path, pattern = paste0(begin, ".*", extension))

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
	
	if (format == "ymd") # year month day
		dt = data.table(file = ls_files, year = ls_files_split[, 1], month = ls_files_split[, 2], day = ls_files_split[, 3])

	if (hour)
	{
		dt[, c("hour", "minute") := as.list(stri_split(str = stri_sub(str = file,
				from = stri_locate_last(file, regex = "_")[, "end"] + 1,
				to = stri_locate_last(file, regex = ".rds")[, "start"] - 1),
			regex = "h", simplify = TRUE)), by = file]
	}

	dt[stri_detect(str = day, regex = extension), day := stri_sub(str = day, to = stri_locate_first(str = day, regex = "_")[,"start"] - 1)]

	setorder(dt, year, month, day, hour, minute)
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
		print(paste0("Figure saved under the name: ", filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf"))
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
lazyPosterior = function(draws, fun = dnorm, filename = NULL, run = NULL, multi = FALSE, ls_nfi = NULL, ...)
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

	# Get the parameter's name if provided
	params = ""
	params_ind = (ls_names == "params")
	if (any(params_ind))
		params = providedArgs[["params"]]

	# Get the index of the x-axis label
	xlab_ind = (ls_names == "xlab")

	# Get the scaling on the x-axis if provided
	scaling_ind = (ls_names == "scaling")
	if (any(scaling_ind)) scaling = providedArgs[["scaling"]] else scaling = 1

	# Get parameters for prior
	if (isTRUE(all.equal(fun, dnorm)))
	{
		if (!all(c("mean", "sd") %in% names(providedArgs)))
			stop("You must provide mean and sd for dnorm")
		
		arg1 = providedArgs[["mean"]]
		arg2 = scaling*providedArgs[["sd"]]
	}

	if (isTRUE(all.equal(fun, dgamma)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape", "rate") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dgamma")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = scaling^2*providedArgs[["var"]] # Squared because in this case it is a variance, not a std. dev.

			arg1 = temp1^2/temp2 # shape
			arg2 = temp1/temp2 # rate
		}

		if (all(c("shape", "rate") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["rate"]]/scaling
		}
	}

	if (isTRUE(all.equal(fun, dbeta)))
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape1", "shape2") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape1 and shape2 for dbeta")

		if (scaling != 1)
			warning("I have not coded the scaling for the beta distribution. Your plot might be out of the window")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = ((1 - temp1)/temp2 - 1/temp1)*temp1^2 # shape 1
			arg2 = arg1*(1/temp1 - 1) # shape 2
		}

		if (all(c("shape1", "shape2") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["shape1"]]
			arg2 = providedArgs[["shape2"]]
		}
		max_y_prior = optimise(f = fun, interval = c(0, 1), maximum = TRUE, shape1 = arg1, shape2 = arg2)[["objective"]]
	}

	# Get posterior
	if (multi)
	{
		info = summary(draws)
		setDT(info)
		length_params = info[, .N]
		density_from_draws = vector(mode = "list", length = length_params)
		x = vector(mode = "list", length = length_params)
		y = vector(mode = "list", length = length_params)
		names(density_from_draws) = info[, variable]
		names(x) = info[, variable]
		names(y) = info[, variable]
		for (varName in info[, variable])
		{
			density_from_draws[[varName]] = density(scaling*draws[, , varName], n = n)
			x[[varName]] = density_from_draws[[varName]]$x
			y[[varName]] = density_from_draws[[varName]]$y
		}
		min_x = min(sapply(x, min))
		max_x = max(sapply(x, max))
		max_y = max(sapply(y, max))
	
	} else {
		density_from_draws = density(draws, n = n)
		x = density_from_draws$x
		y = density_from_draws$y
		min_x = min(x)
		max_x = max(x)
		max_y = max(y)
	}

	min_x = ifelse (min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
	max_x = ifelse (max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x

	if (isTRUE(all.equal(fun, dnorm)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE, mean = arg1, sd = arg2)[["objective"]]
		check_min_bound = integrate(fun, lower = ifelse (min_x < 0, 10*min_x, -10*min_x), upper = min_x, mean = arg1, sd = arg2,
			subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		while (check_min_bound$value > 0.1)
		{
			min_x = ifelse (min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
			check_min_bound = integrate(fun, lower = ifelse (min_x < 0, 10*min_x, -10*min_x), upper = min_x, mean = arg1, sd = arg2,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		}

		check_max_bound = integrate(fun, lower = max_x, upper = ifelse (max_x < 0, -10*max_x, 10*max_x), mean = arg1, sd = arg2,
			subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		while (check_max_bound$value > 0.1)
		{
			max_x = ifelse (max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
			check_max_bound = integrate(fun, lower = max_x, upper = ifelse (max_x < 0, -10*max_x, 10*max_x), mean = arg1, sd = arg2,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		}
	}

	if (isTRUE(all.equal(fun, dgamma)))
	{
		max_y_prior = optimise(f = fun, interval = c(min_x, max_x), maximum = TRUE, shape = arg1, rate = arg2)[["objective"]]
		check_min_bound = integrate(fun, lower = ifelse (min_x < 0, 10*min_x, -10*min_x), upper = min_x, shape = arg1, rate = arg2,
			subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		while (check_min_bound$value > 0.1)
		{
			min_x = ifelse (min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
			check_min_bound = integrate(fun, lower = ifelse (min_x < 0, 10*min_x, -10*min_x), upper = min_x, shape = arg1, rate = arg2,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		}

		check_max_bound = integrate(fun, lower = max_x, upper = ifelse (max_x < 0, -10*max_x, 10*max_x), shape = arg1, rate = arg2,
			subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		while (check_max_bound$value > 0.1)
		{
			max_x = ifelse (max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
			check_max_bound = integrate(fun, lower = max_x, upper = ifelse (max_x < 0, -10*max_x, 10*max_x), shape = arg1, rate = arg2,
				subdivisions = 2000, rel.tol = .Machine$double.eps^0.1)
		}
	}
	max_y = max(max_y, max_y_prior)
	max_y = ifelse (max_y < 0, 0.9*max_y, 1.1*max_y) # To extend 10% from max_y

	# Plot
	if (!is.null(filename))
	{
		filename = paste0(filename, ifelse(!is.null(run), paste0("_", run), ""), ".pdf")
		pdf(filename)
		print(paste0("Figure saved under the name: ", filename))
	}
	
	# Plot posterior
	if (multi)
	{
		colours = MetBrewer::met.brewer("Hokusai3", length_params)
		colours_str = grDevices::colorRampPalette(colours)(length_params)
		colours_str_pol = paste0(colours_str, "66")
		plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", main = paste("Prior and posterior", params),
			xlab = ifelse(any(xlab_ind), providedArgs[["xlab"]], ""))
		for (i in 1:length_params)
		{
			lines(x = density_from_draws[[i]]$x, y = density_from_draws[[i]]$y, col = colours_str[i], lwd = 2)
			polygon(density_from_draws[[i]], col = colours_str_pol[i])
		}
	} else {
		plot(density_from_draws, xlim = c(min_x, max_x), col = "#295384", lwd = 2, main = paste("Prior and posterior", params))
		polygon(density_from_draws, col = "#29538466")
	}

	# Plot prior
	curve(fun(x, arg1, arg2), add = TRUE, lwd = 2, col = "#F4C430")
	DescTools::Shade(fun(x, arg1, arg2), breaks = c(min_x, max_x), col = "#F4C43066", density = NA)

	# Add legend
	if (multi & !is.null(ls_nfi))
		if (length(ls_nfi) != length_params)
			warning("Dimension mismatch between ls_nfi and length_params! The legend might not be correctly printed")

	if (!multi & !is.null(ls_nfi))
		if (length(ls_nfi) != 1)
			warning("To many NFI provided in ls_nfi! The legend might not be correctly printed")

	if (multi)
	{
		legend(x = "topleft", legend = c("Prior", paste("Posterior", if (!is.null(ls_nfi)) ls_nfi else 1:length_params)),
			fill = c("#F4C430", colours_str), box.lwd = 0)
	} else {
		legend(x = "topleft", legend = c("Prior", paste("Posterior", ifelse(!is.null(ls_nfi), ls_nfi, ""))),
			fill = c("#F4C430", "#295384"), box.lwd = 0)
	}

	if (!is.null(filename))
		dev.off()

	return(list(arg1 = arg1, arg2 = arg2, min_x = min_x, max_x = max_x, max_y = max_y, max_y_prior = max_y_prior, filename = filename))
}

## Function to expand the basic names when there is more than one NFI
expand = function(base_names, nb_nfi, patterns = c("Obs", "proba"))
{
	if (nb_nfi < 1)
		stop("Nothing to expand")
	
	new_names = vector(mode = "list", length = length(patterns))
	old_names = base_names
	for (i in 1:length(patterns))
	{
		reg = patterns[i]
		toModify = base_names[stri_detect(base_names, regex = reg)]
		base_names = base_names[!stri_detect(base_names, regex = reg)]
		new_names[[i]] = character(length = nb_nfi*length(toModify))
		for (j in 1:length(toModify))
			new_names[[i]][((j - 1)*nb_nfi + 1):(j*nb_nfi)] = paste0(toModify[j], "[", 1:nb_nfi, "]")
	}

	combined_names = c(base_names, unlist(new_names))
	return(list(old_names = old_names, new_names = combined_names))
}

## Function to do a pair plot of parameters versus energy
energyPairs = function(path, run, results, nb_nfi, energy, rm_names = "latent", n_rows = 3, n_cols = 6, filename = "pairs.pdf")
{
	filename = paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "pairs.pdf")

	params_names = results$metadata()$stan_variables
	params_names = params_names[params_names != "averageGrowth"]
	for (i in 1:length(rm_names))
		params_names = params_names[!stri_detect(params_names, regex = rm_names[i])]
	
	if (nb_nfi > 1)
		params_names = expand(params_names, nb_nfi)[["new_names"]]

	if (length(params_names) > n_rows*n_cols)
		warning("Not enough rows or cols for the number of provided parameters")
	
	pdf(filename, height = 10, width = 10)
	par(mfrow = c(n_rows, n_cols))

	length_params = length(params_names)
	correl_energy = data.table(parameters = params_names, correlation_energy = numeric(length_params))

	for (i in 1:length_params)
	{
		smoothScatter(x = energy, y = as.vector(results$draws(params_names[i])), ylab = params_names[i])
		correl_energy[i, correlation_energy := cor(energy, as.vector(results$draws(params_names[i])))]
	}
	dev.off()

	print(paste("Figure saved under the name:", filename))

	return (correl_energy)
}

## Detect if a species has been processed or not
isProcessed = function(path, multi, lim_time, begin = "^growth-", extension = ".rds$", format = "ymd", lower = 1, upper = 4)
{
	if (class(lim_time) != "Date")
		lim_time = as.Date(lim_time)
	run_vec = lower:upper
	n_runs = length(run_vec)
	processed = rep(FALSE, n_runs)

	if (!dir.exists(path))
		return (FALSE)

	if (multi)
	{
		for (i in 1:n_runs)
		{
			run = run_vec[i]
			date_run = getLastRun(path, begin = begin, extension = extension, format = format, run = i)[["time_ended"]]
			date_run = as.Date(date_run)
			if (!is.na(date_run) & (date_run > lim_time))
				processed[i] = TRUE
		}
	} else {
		date_run = getLastRun(path, begin = begin, extension = extension, format = format, run = 1)[["time_ended"]]
		date_run = as.Date(date_run)
		if (!is.na(date_run) & (date_run > lim_time))
			processed = TRUE
	}
	return (all(processed))
}

## Function to call all the other plot functions and to gather informations on the runs
centralised_fct = function(species, multi, n_runs, ls_nfi, run = NULL, isDBH_normalised = TRUE)
{
	path = paste0("./", species, "/")
	n_nfi = length(ls_nfi)
	error_dt = data.table(nfi = rep(c(ls_nfi, "all")), sigmaProc = numeric(n_nfi + 1), sigmaObs = numeric(n_nfi + 1),
		etaObs = numeric(n_nfi + 1), proba = numeric(n_nfi + 1), correl_sigma_eta = numeric(n_nfi + 1),
		correl_sigma_proba = numeric(n_nfi + 1), correl_eta_proba = numeric(n_nfi + 1))
	setkey(error_dt, nfi)
	
	if (multi & !is.null(run))
	{
		print(paste0("Run = ", run, "; multi and n_runs parameters ignored"))
		temporary = centralised_fct(species, FALSE, n_runs, ls_nfi, run) # Recursive call
		error_dt = temporary[["error_dt"]]
		correl_energy = temporary[["correl_energy"]]
	}

	if (multi & is.null(run))
	{
		error_ls = vector(mode = "list", length = n_runs)
		correl_ls = vector(mode = "list", length = n_runs)
		for (i in 1:n_runs)
		{
			temporary = centralised_fct(species, FALSE, n_runs, ls_nfi, i) # Recursive call
			error_ls[[i]] = temporary[["error_dt"]]
			correl_ls[[i]] = temporary[["correl_energy"]]
		}
		error_dt = rbindlist(error_ls) # , idcol = "run_id"
		correl_energy = rbindlist(correl_ls) # , idcol = "run_id"
	} else {
		# Get inf last run
		info_lastRun = getLastRun(path = path, run = run)
		lastRun = info_lastRun[["file"]]
		time_ended = info_lastRun[["time_ended"]]

		# Load dbh standardising data if necessary
		sd_dbh = 1
		if (isDBH_normalised)
		{
			norm_dbh_dt = readRDS(paste0(path, ifelse(is.null(run), "", paste0(run, "_")), "dbh_normalisation.rds"))
			sd_dbh = norm_dbh_dt[, sd]
		}

		# Load results and associated data set
		results = readRDS(paste0(path, lastRun))
		# stanData = readRDS(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))

		results$print(c("lp__", "averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
			"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), max_rows = 20)

		## Pairs plot of parameters versus energy
		# Extract energy 
		energy = as.vector(results$sampler_diagnostics(inc_warmup = FALSE)[, , "energy__"])
		length(energy)

		# Plot
		correl_energy = energyPairs(path, run = run, results, nb_nfi = n_nfi, energy)

		## Plot prior and posterior for error terms
		# For process error (sigmaProc), it is a variance (of a gamma distrib)
		sigmaProc_array = results$draws("sigmaProc")
		lazyPosterior(draws = sigmaProc_array, fun = dgamma, filename = paste0(path, "sigmaProc_posterior"), params = "process error",
			shape = 5.0^2/1, rate = sd_dbh^2*5.0/1, run = run)

		error_dt["all", sigmaProc := sd_dbh*sqrt(mean(sigmaProc_array))]

		# For measurement error:
		# --- Common variable
		multi_NFI = if (n_nfi > 1) TRUE else FALSE
		# --- Routine measurement error (sigmaObs), it is a sd (of a normal distrib)
		sigmaObs_array = results$draws("sigmaObs")
		lazyPosterior(draws = sigmaObs_array, fun = dgamma, filename = paste0(path, "sigmaObs_posterior"), run = run, xlab = "Error in mm",
			shape = 3.0/0.025, rate = sd_dbh*sqrt(3.0)/0.025, params = "routine obs error", multi = multi_NFI, scaling = sd_dbh,
			ls_nfi = ls_nfi)

		# --- Extreme error (etaObs), it is a sd (of a normal distrib)...
		etaObs_array = results$draws("etaObs")
		lazyPosterior(draws = etaObs_array, fun = dgamma, filename = paste0(path, "etaObs_posterior"), run = run, xlab = "Error in mm",
			shape = 25.6^2/6.2, rate = sd_dbh*25.6/6.2, params = "extreme obs error", multi = multi_NFI, scaling = sd_dbh,
			ls_nfi = ls_nfi)

		# --- ... and its associated probability of occurrence
		proba_array = results$draws("proba")
		lazyPosterior(draws = proba_array, fun = dbeta, filename = paste0(path, "proba_posterior"), run = run, xlab = "Probability",
			shape1 = 48.67, shape2 = 1714.84, params = "probability extreme obs error", multi = multi_NFI,
			ls_nfi = ls_nfi)

		error_dt["all", sigmaObs := sd_dbh*mean(sigmaObs_array)]
		error_dt["all", etaObs := sd_dbh*mean(etaObs_array)]
		error_dt["all", proba := mean(proba_array)]
		error_dt["all", correl_sigma_eta := cor(sigmaObs_array, etaObs_array)]
		error_dt["all", correl_sigma_proba := cor(sigmaObs_array, proba_array)]
		error_dt["all", correl_eta_proba := cor(etaObs_array, proba_array)]

		count = 0

		for (nfi in ls_nfi)
		{
			count = count + 1
			sigmaObs_array = results$draws(paste0("sigmaObs[", count, "]"))
			etaObs_array = results$draws(paste0("etaObs[", count, "]"))
			proba_array = results$draws(paste0("proba[", count, "]"))

			error_dt[nfi, sigmaProc := NA]
			error_dt[nfi, sigmaObs := sd_dbh*mean(sigmaObs_array)]
			error_dt[nfi, etaObs := sd_dbh*mean(etaObs_array)]
			error_dt[nfi, proba := mean(proba_array)]
			error_dt[nfi, correl_sigma_eta := cor(sigmaObs_array, etaObs_array)]
			error_dt[nfi, correl_sigma_proba := cor(sigmaObs_array, proba_array)]
			error_dt[nfi, correl_eta_proba := cor(etaObs_array, proba_array)]
		}

		## Plot prior and posterior for main parameters
		lazyPosterior(draws = results$draws("averageGrowth_mu", inc_warmup = FALSE), fun = dnorm,
			filename = paste0(path, "averageGrowth_mu"), run = run, params = "Average growth (mu)", mean = 0, sd = 1000)
		lazyPosterior(draws = results$draws("averageGrowth_sd", inc_warmup = FALSE), fun = dgamma,
			filename = paste0(path, "averageGrowth_sd"), run = run, params = "Average growth (sd)", shape = 1/100, rate = 1/100)
		
		lazyPosterior(draws = results$draws("dbh_slope", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "dbh_slope"),
			run = run, params = "Dbh slope", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("pr_slope", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "pr_slope"),
			run = run, params = "Precipitation slope", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("pr_slope2", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "pr_slope2"),
			run = run, params = "Precipitation slope (quadratic term)", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("tas_slope", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "tas_slope"),
			run = run, params = "Temperature slope", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("tas_slope2", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "tas_slope2"),
			run = run, params = "Temperature slope (quadratic term)", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("ph_slope", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "ph_slope"),
			run = run, params = "Soil acidity slope (pH)", mean = 0, sd = 5)
		lazyPosterior(draws = results$draws("ph_slope2", inc_warmup = FALSE), fun = dnorm, filename = paste0(path, "ph_slope2"),
			run = run, params = "Soil acidity slope (pH, quadratic term)", mean = 0, sd = 5)

		lazyPosterior(draws = results$draws("competition_slope", inc_warmup = FALSE), fun = dnorm,
			filename = paste0(path, "competition_slope"), run = run, params = "Competition slope", mean = 0, sd = 5)

		## Add run id to data tables
		error_dt[, run_id := run]
		correl_energy[, run_id := run]
	}
	return (list(error_dt = error_dt, correl_energy = correl_energy))
}

## Function to plot the correlations energy <--> parameters for each species, and to plot rescaled error values
plot_correl_error = function(error_dt, correl_dt, threshold_correl = 0.1, rm_correl = "lp__")
{
	correl_dt = correl_dt[!(parameters %in% rm_correl)]
	
	ls_species = correl_dt[, unique(speciesName_sci)] # Should be the same species in error

	for (species in ls_species)
	{
		path = paste0("./", species, "/")
		current_correl = correl_dt[speciesName_sci == species]
		ls_runs = current_correl[speciesName_sci == species, unique(run_id)]
		species_runs = length(ls_runs)
		ls_params = current_correl[speciesName_sci == species, unique(parameters)]
		nb_params = length(ls_params)

		# Plot correlations parameters versus energy
		pdf(paste0(path, "correlations_energy.pdf"), height = 10, width = 10)
		par(mar = c(5, 10, 0, 0))
		plot(0, type = "n", xlim = c(-1, 1), ylim = c(0, nb_params + 1), ylab = "", xlab = "Correlation with energy",
			axes = FALSE, xaxt = "n", yaxt = "n")

		param_counter = 1
		for (current_param in ls_params)
		{
			mean_value = current_correl[parameters == current_param, mean(correlation_energy)]
			# "#72BCD5"
			segments(x0 = 0, y0 = param_counter, x1 = mean_value, y1 = param_counter, lty = 1, lwd = 2,
				col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
			points(mean_value, y = param_counter, pch = 15, cex = 1,
				col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
			run_counter = 1
			if (species_runs %% 2 == 0) # Note that since species_runs > 0, this also implies that species_runs > 1
			{
				intercept_vec = seq(-species_runs/2, species_runs/2, by = 1)
				intercept_vec = intercept_vec[intercept_vec != 0]
				for (i in intercept_vec)
				{
					value = current_correl[parameters == current_param & run_id == run_counter, correlation_energy]
					intercept = param_counter + 0.1*i
					segments(x0 = 0, y0 = intercept, x1 = value, y1 = intercept, lty = 2, lwd = 0.75,
						col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
					points(value, y = intercept, pch = 20, cex = 0.75,
						col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
					run_counter = run_counter + 1
				}
			}
			param_counter = param_counter + 1
		}
		x = c(-threshold_correl, threshold_correl)
		y = c(nb_params + species_runs/2*0.1, nb_params + species_runs/2*0.1)
		lines(rep(x, each = 3), t(matrix(c(y, rep(c(par()$usr[3], NA), each = 2)), ncol = 3)), lwd = 0.5, lty = "dashed")
		lines(rep(0, each = 3), t(matrix(c(nb_params + species_runs/2*0.1, rep(c(par()$usr[3], NA), each = 1)), ncol = 3)), lwd = 2)
		axis(side = 1, at = c(-1, -0.5, -threshold_correl, 0, threshold_correl, 0.5, 1))
		axis(side = 2, at = 1:nb_params, labels = ls_params, las = 1)
		dev.off()

		# Plot correlations among errors and proba
		current_error = error_dt[speciesName_sci == species]
		ls_cols = c("correl_sigma_eta", "correl_sigma_proba", "correl_eta_proba")
		ls_nfi = current_error[nfi != "all", unique(nfi)]
		n_nfi = length(ls_nfi)
		col_counter = 1
		pdf(paste0(path, "correlations_errors_proba.pdf"), height = 10, width = 10)
		if (n_nfi > 1)
		{
			layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2))
			par(mar = c(5, 10, 0, 0))
			plot(0, type = "n", xlim = c(-1, 1), ylim = c(0, length(ls_cols)*n_nfi + 1), ylab = "", xlab = "Correlation",
				axes = FALSE, xaxt = "n", yaxt = "n")
			current_col = ls_cols[1]
			for (current_col in ls_cols)
			{
				for (nfi_id in 1:n_nfi)
				{
					mean_value = mean(unlist(current_error[nfi == ls_nfi[nfi_id], ..current_col]))
					segments(x0 = 0, y0 = col_counter, x1 = mean_value, y1 = col_counter, lty = 1, lwd = 2,
						col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
					points(mean_value, y = col_counter, pch = 15, cex = 1,
						col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
					run_counter = 1
					if (species_runs %% 2 == 0) # Note that this since species_runs > 0, it also implies species_runs > 1 in this case
					{
						intercept_vec = seq(-species_runs/2, species_runs/2, by = 1)
						intercept_vec = intercept_vec[intercept_vec != 0]
						for (i in intercept_vec)
						{
							value = unlist(current_error[nfi == ls_nfi[nfi_id] & run_id == run_counter, ..current_col])
							intercept = col_counter + 0.1*i
							segments(x0 = 0, y0 = intercept, x1 = value, y1 = intercept, lty = 2, lwd = 0.75,
								col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
							points(value, y = intercept, pch = 20, cex = 0.75,
								col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
							run_counter = run_counter + 1
						}
					}
					col_counter = col_counter + 1
				}
			}
			x = c(-threshold_correl, threshold_correl)
			y = c(length(ls_cols)*n_nfi, length(ls_cols)*n_nfi) + species_runs/2*0.1
			lines(rep(x, each = 3), t(matrix(c(y, rep(c(par()$usr[3], NA), each = 2)), ncol = 3)), lwd = 0.5, lty = "dashed")
			lines(rep(0, each = 3), t(matrix(c(length(ls_cols)*n_nfi + species_runs/2*0.1,
				rep(c(par()$usr[3], NA), each = 1)), ncol = 3)), lwd = 2)
			axis(side = 1, at = c(-1, -0.5, -threshold_correl, 0, threshold_correl, 0.5, 1))

			y_labels = expand.grid(c("sigma - eta", "sigma - proba", "eta - proba"), ls_nfi, stringsAsFactors = FALSE)
			setDT(y_labels)
			setorder(y_labels)
			order_according_cols = c((n_nfi + 1):(2*n_nfi), (2*n_nfi + 1):(3*n_nfi), 1:n_nfi)
			axis(side = 2, at = 1:(length(ls_cols)*n_nfi), las = 1,
				labels = do.call(paste, y_labels[order_according_cols]))
		}
		col_counter = 1
		par(mar = c(5, 10, 2, 0))
		plot(0, type = "n", xlim = c(-1, 1), ylim = c(0, length(ls_cols)*n_nfi + 1), ylab = "", xlab = "Correlation",
			axes = FALSE, xaxt = "n", yaxt = "n")
		current_col = ls_cols[1]
		for (current_col in ls_cols)
		{
			mean_value = mean(unlist(current_error[nfi == "all", ..current_col]))
			segments(x0 = 0, y0 = 1.5*col_counter, x1 = mean_value, y1 = 1.5*col_counter, lty = 1, lwd = 2,
				col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
			points(mean_value, y = 1.5*col_counter, pch = 15, cex = 1,
				col = if (abs(mean_value) > threshold_correl) "#EF8A47" else "#34568B")
			run_counter = 1
			if (species_runs %% 2 == 0) # Note that this since species_runs > 0, it also implies species_runs > 1 in this case
			{
				intercept_vec = seq(-species_runs/2, species_runs/2, by = 1)
				intercept_vec = intercept_vec[intercept_vec != 0]
				for (i in intercept_vec)
				{
					value = unlist(current_error[nfi == "all" & run_id == run_counter, ..current_col])
					intercept = 1.5*col_counter + 0.1*i
					segments(x0 = 0, y0 = intercept, x1 = value, y1 = intercept, lty = 2, lwd = 0.75,
						col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
					points(value, y = intercept, pch = 20, cex = 0.75,
						col = if (abs(value) > threshold_correl) "#EF8A47" else "#34568B")
					run_counter = run_counter + 1
				}
			}
			col_counter = col_counter + 1
		}
		x = c(-threshold_correl, threshold_correl)
		y = c(length(ls_cols)*n_nfi, length(ls_cols)*n_nfi) + species_runs/2*0.1
		lines(rep(x, each = 3), t(matrix(c(y, rep(c(par()$usr[3], NA), each = 2)), ncol = 3)), lwd = 0.5, lty = "dashed")
		lines(rep(0, each = 3), t(matrix(c(length(ls_cols)*n_nfi + species_runs/2*0.1,
			rep(c(par()$usr[3], NA), each = 1)), ncol = 3)), lwd = 2)
		axis(side = 1, at = c(-1, -0.5, -threshold_correl, 0, threshold_correl, 0.5, 1))
		axis(side = 2, at = 1.5*1:3, las = 1,
			labels = c("sigma - eta", "sigma - proba", "eta - proba"))

		dev.off()
		print(paste(species, "done"))
	}
}

#### Common variables
## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

n_runs = 4 # Number of runs used in growth_subsample.R
threshold_indiv = 8000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2022/05/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(speciesName_sci, multiRun, threshold_time, lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

error_ls = vector(mode = "list", length = infoSpecies[, .N])
names(error_ls) = infoSpecies[, speciesName_sci]

correl_ls = vector(mode = "list", length = infoSpecies[, .N])
names(correl_ls) = infoSpecies[, speciesName_sci]

#### For loop on processed species, to plot posteriors of the main parameters (errors, intercept, and slopes)
for (species in infoSpecies[, speciesName_sci])
{
	multi = infoSpecies[species, multiRun]
	ls_nfi = unlist(stri_split(infoSpecies[species, ls_nfi], regex = ", "))
	summary_dt = centralised_fct(species, multi, n_runs, ls_nfi, run = if (multi) NULL else 1)
	error_ls[[species]] = summary_dt[["error_dt"]]
	correl_ls[[species]] = summary_dt[["correl_energy"]]
}

error_dt = rbindlist(error_ls, idcol = "speciesName_sci")
correl_dt = rbindlist(correl_ls, idcol = "speciesName_sci")

saveRDS(error_dt, "./error_species.rds")
saveRDS(correl_dt, "./correlation_energy_species.rds")

plot_correl_error(error_dt, correl_dt, threshold_correl = 0.2, rm_correl = "lp__")


#* --------------------------------------------------------------------------------------------
#! ********************************************************************************************
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================    WHAT FOLLOWS ARE OLDER STUFF, SOME ARE GOOD    ================ ODOT#
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! ********************************************************************************************
#* --------------------------------------------------------------------------------------------


#### Plot chains main parameters
lazyTrace(draws = results$draws("averageGrowth_mu", inc_warmup = FALSE), filename = paste0(path, "averageGrowth_mu"), run = run)
lazyTrace(draws = results$draws("averageGrowth_sd", inc_warmup = FALSE), filename = paste0(path, "averageGrowth_sd"), run = run)
lazyTrace(draws = results$draws("dbh_slope", inc_warmup = FALSE), filename = paste0(path, "dbh_slope"), run = run)

lazyTrace(draws = results$draws("pr_slope", inc_warmup = FALSE), filename = paste0(path, "pr_slope"), run = run)
lazyTrace(draws = results$draws("pr_slope2", inc_warmup = FALSE), filename = paste0(path, "pr_slope2"), run = run)
lazyTrace(draws = results$draws("tas_slope", inc_warmup = FALSE), filename = paste0(path, "tas_slope"), run = run)
lazyTrace(draws = results$draws("tas_slope2", inc_warmup = FALSE), filename = paste0(path, "tas_slope2"), run = run)
lazyTrace(draws = results$draws("ph_slope", inc_warmup = FALSE), filename = paste0(path, "ph_slope"), run = run)
lazyTrace(draws = results$draws("ph_slope2", inc_warmup = FALSE), filename = paste0(path, "ph_slope2"), run = run)

lazyTrace(draws = results$draws("competition_slope", inc_warmup = FALSE), filename = paste0(path, "competition_slope"), run = run)

lazyTrace(draws = results$draws("sigmaProc", inc_warmup = FALSE), filename = paste0(path, "sigmaProc"), run = run)

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

# pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "residuals_obs_hist.pdf"))
# hist(dt_dharma[, residuals_obs])
# dev.off()

# qq = qqnorm(dt_dharma[, residuals_obs], pch = 1, frame = FALSE, plot.it = FALSE)

# jpeg(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "qqplot.jpg"))
# plot(qq)
# qqline(dt_dharma[, residuals_obs], col = "#34568B", lwd = 3)
# dev.off()

## Check that the process error does not show any pattern with any predictor
latentG_residuals_array = generate_quantities$draws("latentG_residuals")
mean(latentG_residuals_array)
sd_dbh^2*mean(latentG_residuals_array)

dim(latentG_residuals_array) # iter_sampling * n_chains * n_latentGrowth

# latentG_residuals_vec = as.vector(latentG_residuals_array) # The order is array[,1,1], array[,2,1], ..., array[,n_chain,1], array[,2,1], ...
# length(latentG_residuals_vec)

# pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "latentG_residuals_check_hist.pdf"))
# hist(latentG_residuals_vec)
# dev.off()

latentG_residuals_avg = apply(latentG_residuals_array, 3, mean) # Average latentG_residuals for each latent growth (or dbh)
length(latentG_residuals_avg)

index_notLastMeasure = 1:n_hiddenState
index_notLastMeasure = index_notLastMeasure[!(index_notLastMeasure %in% indices[stanData$last_child_index, index_gen])]

if (length(index_notLastMeasure) != length(latentG_residuals_avg))
	stop("Dimension mismatch between latentG_residuals_avg and index_notLastMeasure")

latent_dbh = apply(generate_quantities$draws("yearly_latent_dbh"), 3, mean)
latent_dbh = latent_dbh[index_notLastMeasure]

pdf(paste0(path, ifelse(!is.null(run), paste0(run, "_"), run), "latentG_residuals_vs_dbh_check.pdf"), height = 6, width = 6)
smoothScatter(x = latent_dbh, y = latentG_residuals_avg)
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
