
#### Aim of prog: Fit the growth data
## Comments:
# 1. To use the GPU with Stan, it is first needed to run clinfo -l. On Bayreuth, I get:
#	clinfo -l
#	Platform #0: NVIDIA CUDA
#		+-- Device #0: NVIDIA RTX A5000
#		`-- Device #1: NVIDIA RTX A5000
#	This indicates that the platform 0 has the GPU with 2 devices.
#	Then, compile the model with the option: cpp_options = list(stan_opencl = TRUE)
#	Finally, to run the model with the function sample, use the argument opencl_ids
#	The opencl_ids is supplied as a vector of length 2, where the first element is the platform ID and the second argument is the device ID.
#	In this case it is (0, 0) or (0, 1) if the second device is desired

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Get parameters for run
args = commandArgs(trailingOnly = TRUE)
# args = c("16", "1", "8000")
if (length(args) != 3)
	stop("Supply the species_id, run_id, and max_indiv as command line arguments!", call. = FALSE)

species_id = as.integer(args[1]) # 17, 48
run_id = as.integer(args[2]) # 1, 2, 3, 4
max_indiv = as.integer(args[3]) # 8000

set.seed(run_id)

#### Tool function
## Initiate Y_gen with reasonable values (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "n_growth", "average_yearlyGrowth", "nbIntervalGrowth", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh_parents, n_growth, average_yearlyGrowth, nbIntervalGrowth, and normalise")

	dbh_parents = providedArgs[["dbh_parents"]]
	n_growth = providedArgs[["n_growth"]]
	average_yearlyGrowth = providedArgs[["average_yearlyGrowth"]]
	nbIntervalGrowth = providedArgs[["nbIntervalGrowth"]]
	normalise = providedArgs[["normalise"]]

	useMean = FALSE

	if ("useMean" %in% names(providedArgs))
		useMean = providedArgs[["useMean"]]

	if (normalise & !all(c("mu_dbh", "sd_dbh") %in% names(providedArgs)))
		stop("You must provide mu_dbh and sd_dbh in order to normalise")

	if (any(average_yearlyGrowth == 0))
	{
		warning("Some average yearly growth were 0. They have been replaced by 0.5")
		average_yearlyGrowth[average_yearlyGrowth == 0] = 0.5
	}

	if (any(average_yearlyGrowth < 0))
	{
		warning("Some average yearly growth were negative. They have been replaced by 0.5")
		average_yearlyGrowth[average_yearlyGrowth < 0] = 0.5
	}

	n_indiv = length(dbh_parents)
	Y_gen = rgamma(n_indiv, dbh_parents^2/0.5, dbh_parents/0.5) # Average = dbh_parents, variance = 0.5

	latent_growth_gen = numeric(n_growth)
	counter_growth = 0

	if (useMean)
	{
		avg_growth = mean(average_yearlyGrowth)
		var_growth = avg_growth/2
		if (avg_growth <= 0)
		{
			warning("Average yearly growth is, in average, negative. Value set to default: 3")
			avg_growth = 3
		}
	}

	# Change extreme growth to more plausible values
	if (!useMean)
	{
		q_90 = quantile(average_yearlyGrowth, seq(0, 1, 0.1))["90%"]
		n_above_q90 = length(average_yearlyGrowth[average_yearlyGrowth > q_90])
		average_yearlyGrowth[average_yearlyGrowth > q_90] = rgamma(n_above_q90, shape = q_90^2/1, rate = q_90/1)

		for (i in 1:n_indiv) # Not that this forbid trees to shrink
		{
			for (j in 1:nbIntervalGrowth[i])
			{
				counter_growth = counter_growth + 1
				latent_growth_gen[counter_growth] = rgamma(n = 1, shape = 2*average_yearlyGrowth[i], rate = 2) # => var = mean/2
			}
		}
	} else {
		for (i in 1:n_indiv) # Not that this forbid trees to shrink
		{
			for (j in 1:nbIntervalGrowth[i])
			{
				counter_growth = counter_growth + 1
				latent_growth_gen[counter_growth] = rgamma(n = 1, shape = avg_growth^2/var_growth, rate = avg_growth/var_growth)
			}
		}
	}

	if (any(latent_growth_gen == 0))
	{
		warning("Some generated latent growth were 0. There have been replaced by 1e-5 (before standardising)")
		latent_growth_gen[latent_growth_gen == 0] = 1e-2
	}

	# Normalise dbh
	if (normalise)
	{
		mu_dbh = providedArgs[["mu_dbh"]]
		sd_dbh = providedArgs[["sd_dbh"]]
		Y_gen = (Y_gen - mu_dbh)/sd_dbh
		latent_growth_gen = latent_growth_gen/sd_dbh
	}

	return(list(latent_dbh_parents = Y_gen, latent_avg_annual_growth = latent_growth_gen))
}

## Function to compute growth, the data table must be sorted by year within tree id and plot id
computeDiametralGrowth = function(dt, col = "growth", byCols = c("plot_id", "tree_id"))
{
	if (!all(c("dbh", "year", byCols) %in% names(dt)))
		stop(paste("The data table must contains at least contains the columns dbh, year,", paste(byCols, collapse = ", ")))
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
			for (i in seq_along(col_start))
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

## Function to subsample the data either spatially (the number of individuals might not be the targeted number) or numerically
subsampling = function(dt, n_indiv_target, mode = "spatial")
{
	if (!any(c("spatial", "numeric") %in% mode))
		stop("Unknown mode, please choose spatial or numeric")

	mean_dbh_beforeSubsample = dt[, mean(dbh)]
	sd_dbh_beforeSubsample = dt[, sd(dbh)]
	quantile_beforeSubsample_25_75 = quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))
	
	if (mode == "spatial")
	{
		n_indiv_per_plot_avg = dt[, length(unique(tree_id)), by = plot_id][, mean(V1)]
		n_plots_sampling = round(n_indiv_target/n_indiv_per_plot_avg)
		coords = dt[, unique(plot_id)]
		sample_plots = sample(x = coords, size = n_plots_sampling)

		dt = dt[plot_id %in% sample_plots]
	}

	if (mode == "numeric")
	{
		parents_index = dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
		sampled_indices = sort(sample(x = parents_index, size = n_indiv_target, replace = FALSE))
		chosen_individuals = dt[sampled_indices, .(plot_id, tree_id)]
		dt = dt[chosen_individuals]
	}

	diffAverage = (dt[, mean(dbh)]/mean_dbh_beforeSubsample > 1.05) | (dt[, mean(dbh)]/mean_dbh_beforeSubsample < 0.95)
	diffSD = (dt[, sd(dbh)]/sd_dbh_beforeSubsample > 1.05) | (dt[, sd(dbh)]/sd_dbh_beforeSubsample < 0.95)
	diffQuantile_25_75 = any((quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))/quantile_beforeSubsample_25_75 > 1.05) |
		(quantile(dt[, dbh], probs = c(0.25, 0.5, 0.75))/quantile_beforeSubsample_25_75 < 0.95))

	n_indiv = unique(dt[, .(tree_id, plot_id)])[, .N]

	return(list(sampledData = dt, diffAverage = diffAverage, diffSD = diffSD, diffQuantile_25_75 = diffQuantile_25_75, n_indiv = n_indiv,
		n_plots_sampling = ifelse(mode == "spatial", n_plots_sampling, NA), mean_dbh_beforeSubsample = mean_dbh_beforeSubsample,
		sd_dbh_beforeSubsample = sd_dbh_beforeSubsample, quantile_beforeSubsample_25_75 = quantile_beforeSubsample_25_75))
}

## Function to average environment
avgClimate = function(climate_dt, indices_avgClim, ls_vars = c("tas", "tasmin", "tasmax", "pr"))
{
	indices_avgClim[, paste0(ls_vars, "_avg") := 0]
	for (i in 1:indices_avgClim[, .N])
		indices_avgClim[i, paste0(ls_vars, "_avg") := climate_dt[index_clim_start:index_clim_end, lapply(.SD, mean), .SDcols = c(ls_vars)]]
}

## Function to recompute indices when subsetting
source("./indices_subsample.R")

#### Load data
## Paths
mainFolder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(mainFolder))
	stop(paste0("Folder\n\t", mainFolder, "\ndoes not exist"))

clim_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(clim_folder))
	stop(paste0("Folder\n\t", clim_folder, "\ndoes not exist"))

soil_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(soil_folder))
	stop(paste0("Folder\n\t", soil_folder, "\ndoes not exist"))

standBasalArea_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(standBasalArea_folder))
	stop(paste0("Folder\n\t", standBasalArea_folder, "\ndoes not exist"))

## Tree inventories data
treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))
setkey(treeData, plot_id, tree_id, year)
ls_species = sort(treeData[, unique(speciesName_sci)])

if ((species_id < 1) | (species_id > length(ls_species)))
	stop(paste0("Species id = ", species_id, " has no corresponding species (i.e., either negative or larger than the number of species)"))

speciesCountry = treeData[, .N, by = .(speciesName_sci, country)]
setkey(speciesCountry, speciesName_sci)
setnames(speciesCountry, "N", "n_measurements")
speciesCountry[, prop := round(100*n_measurements/sum(n_measurements), 2), by = speciesName_sci]

species = ls_species[species_id]
print(paste("Script running for species:", species))
print(speciesCountry[species])

treeData = treeData[speciesName_sci == species]

savingPath = paste0("./", species, "/")

if (!dir.exists(savingPath))
	dir.create(savingPath)

## Subsample tree data if necessary
n_indiv = unique(treeData[, .(tree_id, plot_id)])[, .N]
print(paste("Number of individuals:", n_indiv))
subsamplingActivated = FALSE

if (n_indiv > max_indiv)
{
	print("Too many individuals, subsampling")
	subsamplingActivated = TRUE

	checkSampling = subsampling(treeData, n_indiv_target = max_indiv, mode = "spatial")
	treeData = checkSampling[["sampledData"]]
	n_indiv = checkSampling[["n_indiv"]]

	if (checkSampling[["diffAverage"]])
		stop("The subsample does not look representative of the whole data set, check the average")

	if (checkSampling[["diffSD"]])
		stop("The subsample does not look representative of the whole data set, check the std. dev")

	if (checkSampling[["diffQuantile_25_75"]])
		stop("The subsample does not look representative of the whole data set, check the quantiles 0.25, 0.5, and 0.75")
}

if ((!subsamplingActivated) & (run_id != 1))
	stop("Running the model only once (i.e., with run_id = 1) is enough: There is no subsampling")

n_inventories = length(treeData[, unique(nfi_id)])

## Compute if there are two measurements or more
if (treeData[, .N, by = .(plot_id, tree_id)][, min(N) < 2])
	stop("There are individuals measured only once")

## Compute growth
computeDiametralGrowth(treeData, byCols = c("speciesName_sci", "plot_id", "tree_id"))
growth_dt = na.omit(treeData)

print(paste0(round(100*growth_dt[growth < 0 | growth > 10, .N]/growth_dt[, .N], 3), "% negative or above 10 mm/yr"))

## Print number of measurements per country
print("Number of growth measurements per country:")
countryStats = growth_dt[, .N, by = country]
countryStats[, prop := round(100*N/sum(N), 2)]
print(countryStats)

## Read climate
climate = readRDS(paste0(clim_folder, "europe_reshaped_climate.rds"))
setkey(climate, plot_id, year)
climate[, row_id := 1:.N]

## Read soil data (pH)
soil = readRDS(paste0(soil_folder, "europe_reshaped_soil.rds"))
soil = soil[plot_id %in% treeData[, plot_id]]

## Read interpolated basal area data
standBasalArea = readRDS(paste0(standBasalArea_folder, "europe_reshaped_standBasalArea.rds"))
setkey(standBasalArea, plot_id, year)

## Set-up indices
indices_list = indices_subsample(run_id, treeData, climate, savingPath, mainFolder, clim_folder)
indices = indices_list[["indices"]]
indices_avgClim = indices_list[["indices_avg"]]

if (indices[, .N] != treeData[, .N])
	stop(paste0("Dimension mismatch between indices and treeData for species `", species, "`"))

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

# Compute the number of measured averaged growth (it is also the number of latent growth)
n_growth = n_obs - n_indiv
nbIntervalGrowth = unique(indices[, .(plot_id, tree_id, nbIntervalGrowth)])[, nbIntervalGrowth]

if (n_growth != sum(nbIntervalGrowth))
	stop("Dimensions mismatch between latent growth and nb intervals growth")

print(paste("Number of latent growth:", n_growth))

# Define parents, children, and last child
parents_index = treeData[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

# Define for each NFI at which individual they start and end (given treeData is sorted by plot_id, with the country first)!
start_nfi_avg_growth = integer(n_inventories)
end_nfi_avg_growth = integer(n_inventories)

ls_countries = treeData[, unique(country)]
start_nfi_avg_growth[1] = 1

if (n_inventories > 1)
{
	for (k in 1:(n_inventories - 1))
	{
		end_nfi_avg_growth[k] = start_nfi_avg_growth[k] + growth_dt[(stri_detect_regex(plot_id, ls_countries[k])), .N] - 1
		start_nfi_avg_growth[k + 1] = end_nfi_avg_growth[k] + 1
	}
}

end_nfi_avg_growth[n_inventories] = n_obs - n_indiv # Which is n_children

if (growth_dt[, .N] != n_obs - n_indiv)
	stop("Dimension mismatch between growth_dt and number of observations")

if (n_inventories == 1)
{
	start_nfi_avg_growth = as.array(start_nfi_avg_growth)
	end_nfi_avg_growth = as.array(end_nfi_avg_growth)
}

## Average climate
avgClimate(climate, indices_avgClim, ls_vars = c("tas", "tasmin", "tasmax", "pr"))
avgClimate(standBasalArea, indices_avgClim, ls_vars = c("standBasalArea_interp"))

#### Compute and save normalising constants
normalisation(dt = treeData, colnames = "dbh", folder = savingPath, filename = paste0(run_id, "_dbh_normalisation_classic.rds"))
normalisation(dt = indices_avgClim, colnames = c("pr_avg", "tas_avg"), folder = savingPath,
	filename = paste0(run_id, "_climate_normalisation_classic.rds"))
normalisation(dt = soil, colnames = "ph", folder = savingPath, filename = paste0(run_id, "_ph_normalisation_classic.rds"))
normalisation(dt = indices_avgClim, colnames = "standBasalArea_interp_avg", folder = savingPath,
	filename = paste0(run_id, "_ba_normalisation_classic.rds"))

climate_mu_sd = readRDS(paste0(savingPath, run_id, "_climate_normalisation_classic.rds"))
ph_mu_sd = readRDS(paste0(savingPath, run_id, "_ph_normalisation_classic.rds"))
ba_mu_sd = readRDS(paste0(savingPath, run_id, "_ba_normalisation_classic.rds"))

#### Stan model
## Define stan variables
# Common variables
maxIter = 2500
n_chains = 3

# Initial values for states only
average_yearlyGrowth = (treeData[last_child_index, dbh] - treeData[parents_index, dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_yearlyGrowth) != n_indiv)
	stop("Dimensions mismatch between average_yearlyGrowth and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, dbh],
	n_growth = n_growth, average_yearlyGrowth = average_yearlyGrowth, nbIntervalGrowth = nbIntervalGrowth,
	normalise = TRUE, mu_dbh = 0, sd_dbh = sd(treeData[, dbh]), useMean = FALSE)

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]][["latent_avg_annual_growth"]]))

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Total number of individuals (all NFIs together)
	n_climate = indices_avgClim[, .N], # Dimension of the climate vector (all NFIs together)
	n_plots = length(treeData[, unique(plot_id)]), # Number of plots (all NFIs together)
	n_obs = n_obs, # Total number of tree observations (all NFIs together)
	n_growth = n_growth, # Number of measured growth = n_obs - n_indiv
	nbIntervalGrowth = nbIntervalGrowth, # Number of growth intervals
	deltaYear = growth_dt[, deltaYear], # Number of years between two measurements of an individual
	n_inventories = n_inventories, # Number of forest inventories involving different measurement errors in the data

	# Indices
	climate_index = unique(indices[, .(tree_id, plot_id, plot_index)])[, plot_index], # Index of the climate associated to each parent
	
	start_nfi_avg_growth = start_nfi_avg_growth, # Starting point of each NFI for averaged obs growth
	end_nfi_avg_growth = end_nfi_avg_growth, # Ending point of each NFI for averaged obs growth
	
	plot_index = unique(indices[, .(tree_id, plot_id, plot_index)])[, plot_index], # Indicates to which plot individuals belong to

	# Observations
	avg_annual_growth_obs = growth_dt[, growth],
	dbh_init = treeData[parents_index, dbh], # Initial dbh data

	# Explanatory variables
	sd_dbh = ifelse(subsamplingActivated, checkSampling[["sd_dbh_beforeSubsample"]], treeData[, sd(dbh)]),

	precip = indices_avgClim[, pr_avg], # Precipitations
	pr_mu = climate_mu_sd[variable == "pr_avg", mu], # To centre the precipitations
	pr_sd = climate_mu_sd[variable == "pr_avg", sd], # To standardise the precipitations

	tas = indices_avgClim[, tas_avg], # Temperature
	tas_mu = climate_mu_sd[variable == "tas_avg", mu], # To centre the temperature
	tas_sd = climate_mu_sd[variable == "tas_avg", sd], # To standardise the temperature

	ph = soil[, ph], # pH of the soil measured with CaCl2
	ph_mu = ph_mu_sd[variable == "ph", mu], # To centre the pH
	ph_sd = ph_mu_sd[variable == "ph", sd], # To standardise the pH

	standBasalArea = indices_avgClim[, standBasalArea_interp_avg], # Average stand basal area for a given plot at a given time interval
	ba_mu = ba_mu_sd[variable == "standBasalArea_interp_avg", mu], # To centre the basal area
	ba_sd = ba_mu_sd[variable == "standBasalArea_interp_avg", sd] # To standardise the basal area
)

saveRDS(object = stanData, file = paste0(savingPath, run_id, "_stanData_classic.rds"))
saveRDS(object = treeData, file = paste0(savingPath, run_id, "_treeData_classic.rds"))

## Compile model
model = cmdstan_model("./growth_classic.stan", cpp_options = list(stan_opencl = TRUE))

start_time = Sys.time()

## Run model
results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 50, chains = n_chains,
	iter_warmup = 1500, iter_sampling = 1000, save_warmup = TRUE, init = initVal_Y_gen,
	max_treedepth = 13, adapt_delta = 0.95, opencl_ids = c(0, ifelse(species_id %% 2 == 0, 0, 1)))

end_time = Sys.time()

time_ended = format(Sys.time(), "%Y-%m-%d_%Hh%M")
results$save_output_files(dir = savingPath, basename = paste0("growth-run=", run_id, "-", time_ended), timestamp = FALSE, random = TRUE)
results$save_object(file = paste0(savingPath, "growth-run=", run_id, "-", time_ended, "_de-fr-sw_", max_indiv, "_main_classic.rds"))

results$cmdstan_diagnose()

print(end_time - start_time)
results$print(c("lp__", "averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "etaObs", "proba", "sigmaProc"), max_rows = 20)
