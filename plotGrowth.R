
#### Aim of prog: Plot different visualisations of growth versus climate and competition

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
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
# shapefile_path = "/Users/mistral/Nextcloud/shapefiles/Germany_shapefile/"

## Tree data
species = "Picea_abies"
path = paste0("./", species, "/")
run = 1
info_lastRun = getLastRun(path = path, run = run)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]
results = readRDS(paste0(path, lastRun))

## Associated estimated parameters
params_names = c("potentialGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "competition_slope",
	"measureError", "processError")

# Mean estimation
estimated_params = getParams(model_cmdstan = results, params_names = params_names)

# Array to create a set of parameters
estimated_params_array = results$draws(params_names)
dim(estimated_params_array) # iter_sampling * n_chain * nbParams

n_rep = results$metadata()$iter_sampling * results$num_chains()
dt_params_set = data.table()
dt_params_set[, c(params_names) := numeric(n_rep)]
for (param in params_names)
	dt_params_set[, (param) := as.vector(estimated_params_array[, , param])]

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
dbh_mu_sd = readRDS(paste0(tree_path, "dbh_normalisation.rds"))
climate_mu_sd = readRDS(paste0(tree_path, "climate_normalisation.rds"))

#### Compute average growth in landscape
average_growth = setDT(as.data.frame(x = climate_mean, xy = TRUE))

lower_dbh = 40
upper_dbh = 700

average_growth[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = precipitation, temperature = temperature, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = .(x, y)]

average_growth[, integral := integral/(upper_dbh - lower_dbh)]

#### Compute average growth along a one-dimension environmental gradient
## Common variables
mean_tas = mean(values(subset(climate_mean, "temperature")), na.rm = TRUE)
min_tas = min(values(subset(climate_mean, "temperature")), na.rm = TRUE)
max_tas = max(values(subset(climate_mean, "temperature")), na.rm = TRUE)

mean_pr = mean(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
min_pr = min(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
max_pr = 1800 # max(values(subset(climate_mean, "precipitation")), na.rm = TRUE)

n_tas = 200
n_pr = 2000

## Temperature gradient
dt_growth_tas = data.table(temperature = seq(min_tas, max_tas, length.out = n_tas), integral = numeric(n_tas))

dt_growth_tas[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = mean_pr, temperature = temperature, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = temperature]

selected_iterations = seq(1, dt_params_set[, .N], by = 3)

for (i in selected_iterations)
{
	current_params = unlist(dt_params_set[i, ])
	dt_growth_tas[, (paste0("integral", i)) := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
		precipitation = mean_pr, temperature = temperature, basalArea = 0,
		params = current_params, scaled_dbh = FALSE, scaled_clim = FALSE,
		sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
		tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = temperature]

	if (i %% 10 == 0)
		print(paste0(i*100/n_rep, "% done"))
}

dt_growth_tas[, integral := integral/(upper_dbh - lower_dbh)]
dt_growth_tas[, paste0("integral", selected_iterations) := lapply(.SD, function(x){x/(upper_dbh - lower_dbh)}),
	.SDcols = paste0("integral", selected_iterations)]

## Precipitation gradient
dt_growth_pr = data.table(precipitation = seq(min_pr, max_pr, length.out = n_pr), integral = numeric(n_pr))

dt_growth_pr[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	precipitation = precipitation, temperature = mean_tas, basalArea = 0,
	params = estimated_params, scaled_dbh = FALSE, scaled_clim = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
	tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = precipitation]

for (i in selected_iterations)
{
	current_params = unlist(dt_params_set[i, ])
	dt_growth_pr[, (paste0("integral", i)) := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
		precipitation = precipitation, temperature = mean_tas, basalArea = 0,
		params = current_params, scaled_dbh = FALSE, scaled_clim = FALSE,
		sd_dbh = dbh_mu_sd[1, sd], pr_mu = climate_mu_sd[variable == "pr", mu], pr_sd = climate_mu_sd[variable == "pr", sd],
		tas_mu = climate_mu_sd[variable == "tas", mu], tas_sd = climate_mu_sd[variable == "tas", sd])$value, by = precipitation]

	if (i %% 10 == 0)
		print(paste0(i*100/n_rep, "% done"))
}

dt_growth_pr[, integral := integral/(upper_dbh - lower_dbh)]
dt_growth_pr[, paste0("integral", selected_iterations) := lapply(.SD, function(x){x/(upper_dbh - lower_dbh)}),
	.SDcols = paste0("integral", selected_iterations)]

#### Plot average growth
## In function of climate
# Common variables
y_min = min(dt_growth_tas[, min(integral)], dt_growth_pr[, min(integral)])
y_max = max(dt_growth_tas[, max(integral)], dt_growth_pr[, max(integral)])

# Average growth versus temperature
# pdf(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "growth_vs_tas.pdf"))
png(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "growth_vs_tas.png"), width = 1080, height = 1080)
par(mar = c(5, 5, 0.5, 0.5))
plot(dt_growth_tas[, temperature], dt_growth_tas[, integral], xlab = "Temperature", ylab = "Growth (mm/yr)",
	type = "l", las = 1, cex.lab = 1.75, cex.axis = 1.75, ylim = c(y_min, y_max), col = "#CD212A", lwd = 4)

for (i in selected_iterations)
{
	col_i = paste0(paste0("integral", i))
	lines(dt_growth_tas[, temperature], unlist(dt_growth_tas[, ..col_i]), col = "#44444422")
}
dev.off()

# Average growth versus precipitation
pdf(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "growth_vs_pr.pdf"))
par(mar = c(5, 5, 0.5, 0.5))
plot(dt_growth_pr[, precipitation], dt_growth_pr[, integral], xlab = "Precipitation", ylab = "Growth (mm/yr)",
	type = "l", las = 1, cex.lab = 1.75, cex.axis = 1.75, ylim = c(y_min, y_max), col = "#CD212A", lwd = 4)

for (i in selected_iterations)
{
	col_i = paste0(paste0("integral", i))
	lines(dt_growth_pr[, precipitation], unlist(dt_growth_precipitation[, ..col_i]), col = "#44444422")
}
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
