
#### Aim of prog: Plot different visualisations of growth versus climate and competition

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(posterior)
library(cmdstanr)
library(stringi)
library(viridis)
library(terra)

#### Tool functions
## Growth function
growth_fct = function(x, precipitation, temperature, basalArea, params, scaled_dbh = FALSE, scaled_clim = FALSE, ...)
{
	providedArgs = list(...)
	providedArgs_names = names(providedArgs)
	scaling = 1

	if (!scaled_dbh)
	{
		if (!("sd_dbh" %in% providedArgs_names))
			stop("You need to provide sd_dbh in order to standardise dbh")
		
		sd_dbh = providedArgs[["sd_dbh"]]
		scaling = sd_dbh
		x = x/sd_dbh
	}

	if (!scaled_clim)
	{
		if (!all(c("pr_mu", "pr_sd", "tas_mu", "tas_sd") %in% providedArgs_names))
			stop("You need to provide pr_mu, pr_sd, tas_mu, and tas_sd in order to standardise dbh")
		
		pr_mu = providedArgs[["pr_mu"]]
		pr_sd = providedArgs[["pr_sd"]]
		tas_mu = providedArgs[["tas_mu"]]
		tas_sd = providedArgs[["tas_sd"]]
		temperature = (temperature - tas_mu)/tas_sd
		precipitation = (precipitation - pr_mu)/pr_sd
	}

	potentialGrowth = params["potentialGrowth"]
	dbh_slope = params["dbh_slope"]

	pr_slope = params["pr_slope"]
	pr_slope2 = params["pr_slope2"]
	tas_slope = params["tas_slope"]
	tas_slope2 = params["tas_slope2"]

	competition_slope = params["competition_slope"]

	return(scaling*exp(potentialGrowth + dbh_slope*x + pr_slope*precipitation + pr_slope2*precipitation^2 +
		tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*basalArea))
}

## Function to get estimated parameters from results (fixed values, either mean or median)
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
		begin = paste0(begin, "run=", run, "-")
	
	ls_files = list.files(path = path, pattern = paste0("^", begin, ".*", extension, "$"))
	ls_files_split = stri_split(stri_sub(ls_files, from = stri_locate(ls_files, regex = begin)[, "end"] + 1),
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

#### Load data
## Common variables
species = "Picea_abies"
run = 1
years = 2010:2018

## Paths
tree_path = paste0("./", species, "/")
climate_path = "/home/amael/project_ssm/climateData/Chelsa/yearlyAverage/"
shapefile_path = "/home/amael/shapefiles/deutschland/"

## Tree data
info_lastRun = getLastRun(path = tree_path, run = run)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]
results = readRDS(paste0(tree_path, lastRun))

## Associated estimated parameters
params_names = c("potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "competition_slope",
	"measureError", "processError")

# Mean estimation
estimated_params = getParams(model_cmdstan = results, params_names = params_names)

## Climate
temperature_files = paste0(climate_path, "tas/", years, ".tif")
precipitation_files = paste0(climate_path, "pr/", years, ".tif")

climate = rast(c(temperature_files, precipitation_files))
names(climate) = c(paste0("tas_", years), paste0("pr_", years))

## Landscape
germany = vect(paste0(shapefile_path, "germany.shp"))
germany = simplifyGeom(x = germany, tolerance = 100)
germany = project(x = germany, y = climate)

bayern = vect(paste0(shapefile_path, "bayern.shp"))
bayern = project(x = bayern, y = climate)

cities = data.table(
	x = c(11.576124, 13.404954, 8.682127, 11.061859),
	y = c(48.137154, 52.520008, 50.110924, 49.460983),
	name = c("München", "Berlin", "Frankfurt", "Nürnberg"))

cities_rs = vect(x = cities, geom = c("x", "y"), crs = crs(germany))

## Crop climate to Germany only
climate = crop(x = climate, y = germany, mask = TRUE)
climate_mean = c(mean(subset(climate, paste0("tas_", years))), mean(subset(climate, paste0("pr_", years))))
names(climate_mean) = c("temperature", "precipitation")

## Scaling
dbh_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "dbh_normalisation.rds"))
climate_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "climate_normalisation.rds"))

#### Compute average growth in a landscape/along a gradient
## Define stan model to compute average growth (integral in generated quantities block)
gq_script = write_stan_file(
"
functions {
	// Function to integrate (growth function). This returns the expected growth in mm, for 1 year.
	real growth_integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i)
	{
		// theta is the vector of parameters
		real potentialGrowth = theta[1];
		real dbh_slope = theta[2];

		real pr_slope = theta[3];
		real pr_slope2 = theta[4];
		real tas_slope = theta[5];
		real tas_slope2 = theta[6];

		real competition_slope = theta[7];

		// x_r is an array of real data (here, precip, competition, etc...)
		real precip = x_r[1];
		real temperature = x_r[2];
		real totalTreeWeight = x_r[3];

		return (exp(potentialGrowth + dbh_slope*x + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_climate_new; // Dimension of the new climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	// array [n_climate_new] real precip_new; // New precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	// array [n_climate_new] real tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0>[n_indiv] totalTreeWeight; // Sum of the tree weights for a given plot at a given time
	// array [n_indiv] real totalTreeWeight_new; // Sum of the tree weights for a given plot at a given time

	array [3*n_climate_new] real x_r; // Contains in this order: pr, ts, totalTreeWeight, each of size n_climate_new, and already standardised
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_indiv] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);

	// Array concatenating the three real predictors (used to compute integral), and unused, but necessary, int array to compute integral
/*	array [3 * n_climate] real x_r;
	for (i in 1:n_climate)
	{
		x_r[i] = normalised_precip[i];
		x_r[i + n_climate] = normalised_tas[i];
		x_r[i + 2*n_climate] = normalised_totalTreeWeight[i];
	}
*/
	array [0] int x_i;
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd(Yobs)^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	array [n_climate_new] real average_growth_sim;
	
	{
		real lower_bound = 50/sd(Yobs);
		real upper_bound = 700/sd(Yobs);
		array [3] int index;

		for (i in 1:n_climate_new)
		{
			index = {i, i + n_climate_new, i + 2*n_climate_new}; // To get the i^th value of precip, temp, and tree weight, respectively
			average_growth_sim[i] = integrate_1d(growth_integrand, lower_bound, upper_bound,
				{potentialGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope},
				x_r[index], x_i)/(upper_bound - lower_bound);
		}
	}
}
"
)

gq_model = cmdstan_model(gq_script)

## Common variables
mean_tas = mean(values(subset(climate_mean, "temperature")), na.rm = TRUE)
min_tas = min(values(subset(climate_mean, "temperature")), na.rm = TRUE)
max_tas = max(values(subset(climate_mean, "temperature")), na.rm = TRUE)

mean_pr = mean(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
min_pr = min(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
max_pr = max(values(subset(climate_mean, "precipitation")), na.rm = TRUE)

n_tas = 200
n_pr = 600

n_chains = results$num_chains()

lower_dbh = 40
upper_dbh = 700

## Compute average growth in a landscape
average_growth = setDT(as.data.frame(x = climate_mean, xy = TRUE))

average_growth[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = precipitation, temperature = temperature, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = .(x, y)]

average_growth[, integral := integral/(upper_dbh - lower_dbh)]

## Compute average growth along a one-dimension environmental gradient
# Temperature gradient
dt_growth_tas = data.table(temperature = seq(min_tas, max_tas, length.out = n_tas), integral = numeric(n_tas))

dt_growth_tas[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = mean_pr, temperature = temperature, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = temperature]

dt_growth_tas[, integral := integral/(upper_dbh - lower_dbh)]

stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
stanData$x_r = c((rep(mean_pr, n_tas) - stanData$pr_mu)/stanData$pr_sd, # Precip
	(dt_growth_tas[, temperature] - stanData$tas_mu)/stanData$tas_sd, # Temperature
	rep(0, n_tas) # totalTreeWeight, note that I set it to 0, which is not the averaged on the centred scale
)
stanData$n_climate_new = n_tas

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

average_growth_sim = generate_quantities$draws("average_growth_sim")
quantile_tas = data.table(q5 = numeric(n_tas), q95 = numeric(n_tas))
for (i in 1:n_tas)
	quantile_tas[i, c("q5", "q95") := as.list(sd(stanData$Yobs)*quantile2(x = average_growth_sim[, , i], probs = c(0.05, 0.95)))]

# Precipitation gradient
dt_growth_pr = data.table(precipitation = seq(min_pr, max_pr, length.out = n_pr), integral = numeric(n_pr))

dt_growth_pr[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = precipitation, temperature = mean_tas, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = precipitation]

dt_growth_pr[, integral := integral/(upper_dbh - lower_dbh)]

# Generate simulations
stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
stanData$x_r = c((dt_growth_pr[, precipitation] - stanData$pr_mu)/stanData$pr_sd, # Precip
	(rep(mean_tas, n_pr) - stanData$tas_mu)/stanData$tas_sd, # Temperature
	rep(0, n_pr) # totalTreeWeight, note that I set it to 0, which is not the averaged on the centred scale
)
stanData$n_climate_new = n_pr

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

average_growth_sim = generate_quantities$draws("average_growth_sim")
quantile_pr = data.table(q5 = numeric(n_pr), q95 = numeric(n_pr))
for (i in 1:n_pr)
	quantile_pr[i, c("q5", "q95") := as.list(sd(stanData$Yobs)*quantile2(x = average_growth_sim[, , i], probs = c(0.05, 0.95)))]

#### Plot average growth
## In function of climate
# Common variables
y_min = min(dt_growth_tas[, min(integral)], dt_growth_pr[, integral],
	quantile_tas[, min(q5)], quantile_pr[, min(q5)])

y_max = max(dt_growth_tas[, max(integral)], dt_growth_pr[, integral],
	quantile_tas[, max(q95)], quantile_pr[, max(q95)])

# Average growth versus temperature
pdf(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "growth_vs_tas.pdf"))
par(mar = c(5, 5, 0.5, 0.5))
plot(dt_growth_tas[, temperature], dt_growth_tas[, integral], xlab = "Temperature", ylab = "Growth (mm/yr)",
	type = "n", las = 1, cex.lab = 1.75, cex.axis = 1.75, ylim = c(y_min, y_max))

polygon(c(rev(dt_growth_tas[, temperature]), dt_growth_tas[, temperature]),
	c(rev(quantile_tas[, q95]), quantile_tas[, q5]), col = '#44444422', border = NA)

lines(dt_growth_tas[, temperature], dt_growth_tas[, integral], col = "#34568B", lwd = 4)

lines(dt_growth_tas[, temperature], quantile_tas[, q95], lty = 'dashed', col = "#F4C430", lwd = 2)
lines(dt_growth_tas[, temperature], quantile_tas[, q5], lty = 'dashed', col = "#F4C430", lwd = 2)

dev.off()

# Average growth versus precipitation
pdf(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "growth_vs_pr.pdf"))
par(mar = c(5, 5, 0.5, 0.5))
plot(dt_growth_pr[, precipitation], dt_growth_pr[, integral], xlab = "Precipitation", ylab = "Growth (mm/yr)",
	type = "n", las = 1, cex.lab = 1.75, cex.axis = 1.75, ylim = c(y_min, y_max))

polygon(c(rev(dt_growth_pr[, precipitation]), dt_growth_pr[, precipitation]),
	c(rev(quantile_pr[, q95]), quantile_pr[, q5]), col = '#44444422', border = NA)

lines(dt_growth_pr[, precipitation], dt_growth_pr[, integral], col = "#34568B", lwd = 4)

lines(dt_growth_pr[, precipitation], quantile_pr[, q95], lty = 'dashed', col = "#F4C430", lwd = 2)
lines(dt_growth_pr[, precipitation], quantile_pr[, q5], lty = 'dashed', col = "#F4C430", lwd = 2)

dev.off()

## In a landscape
growth_rs = rast(x = average_growth, type = "xyz", crs = crs(climate_mean))

png(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "average_growth_landscape.png"), width = 1080, height = 1080)
plot(subset(growth_rs, "integral"), axes = FALSE, col = viridis(256),
	plg = list(cex = 3, title = "Growth(mm/yr)"))
plot(germany, add = TRUE, col = NA, lwd = 2)
plot(cities_rs, pch = 19, add = TRUE, cex = 4)
dev.off()

## In Bavaria
growth_bayern_rs = crop(growth_rs, bayern, mask = TRUE)

png(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "average_growth_landscape_Bayern.png"), width = 1080, height = 1080)
plot(subset(growth_bayern_rs, "integral"), axes = FALSE, col = viridis(256),
	plg = list(cex = 3, title = "Growth(mm/yr)"))
plot(bayern, add = TRUE, lwd = 2)
plot(cities_rs, pch = 19, add = TRUE, cex = 4)
dev.off()
