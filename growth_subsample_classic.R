
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
	requiredArgs = c("dbh_parents", "n_latentGrowth", "average_yearlyGrowth", "nbIntervalGrowth", "normalise")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide dbh_parents, n_latentGrowth, average_yearlyGrowth, nbIntervalGrowth, and normalise")

	dbh_parents = providedArgs[["dbh_parents"]]
	n_latentGrowth = providedArgs[["n_latentGrowth"]]
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

	latent_growth_gen = numeric(n_latentGrowth)
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
		coords = unique(dt[, plot_id])
		sample_plots = sample(x = coords, size = n_plots_sampling)

		setkey(dt, plot_id)
		dt = dt[sample_plots]
		setorder(dt, plot_id, tree_id, year)
	}

	if (mode == "numeric")
	{
		setkey(dt, tree_id, plot_id)
		parents_index = dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
		sampled_indices = sort(sample(x = parents_index, size = n_indiv_target, replace = FALSE))
		chosen_individuals = dt[sampled_indices, .(tree_id, plot_id)]
		dt = dt[chosen_individuals]
		setorder(dt, plot_id, tree_id, year)
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

## Function to recompute indices when subsetting
source("./indices_subsample_classic.R")

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

## Read soil data (pH)
soil = readRDS(paste0(soil_folder, "europe_reshaped_soil.rds"))

## Read interpolated basal area data
standBasalArea = readRDS(paste0(standBasalArea_folder, "europe_reshaped_standBasalArea.rds"))

## Average climate

