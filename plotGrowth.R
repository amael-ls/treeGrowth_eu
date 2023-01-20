
#### Aim of prog: Plot different visualisations of growth versus predictors (diameters, environmental factors)

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
source("toolFunctions.R")

#### Load data
## Common variables
species = "Abies alba"
run = 1
years = 2010:2018

## Paths
tree_path = paste0("./", species, "/")
climate_path = "/bigdata/Predictors/climateChelsa/"
ph_path = "/bigdata/Predictors/Soil\ esdacph\ Soil\ pH\ in\ Europe/"
shapefile_path = "/home/amael/shapefiles/Deutschland/"

## Results
# State-Space Model approach
info_lastRun = getLastRun(path = tree_path, begin = "growth-", extension = "_main.rds$", run = run)
(lastRun = info_lastRun[["file"]])
results = readRDS(paste0(tree_path, lastRun))
nb_nfi = results$metadata()$stan_variable_sizes$etaObs

# Classic approach
info_lastRun = getLastRun(path = tree_path, begin = "growth-", extension = "_classic.rds$", run = run)
(lastRun = info_lastRun[["file"]])
results_classic = readRDS(paste0(tree_path, lastRun))

## Associated estimated parameters
params_names = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "etaObs", "proba", "sigmaProc")

if (nb_nfi > 1)
	params_names = expand(params_names, nb_nfi)[["new_names"]]

# Mean estimation
estimated_params = getParams(model_cmdstan = results, params_names = params_names)

## Climate (now, it would be more accurate to call this predictors, but by that time there were only climate!)
temperature_files = paste0(climate_path, "tas/", years, ".tif")
precipitation_files = paste0(climate_path, "pr/", years, ".tif")

soil_shp = vect(paste0(ph_path, "country_laea.shp"))
soil = rast(paste0(ph_path, "ph_cacl2/w001001.adf"))
crs(soil) = crs(soil_shp)

climate = rast(c(temperature_files, precipitation_files))
soil = project(soil, climate)
climate = c(climate, soil)
names(climate) = c(paste0("tas_", years), paste0("pr_", years), "ph")

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
climate_mean = c(mean(subset(climate, paste0("tas_", years))), mean(subset(climate, paste0("pr_", years))), mean(subset(climate, "ph")))
names(climate_mean) = c("temperature", "precipitation", "ph")

## Scaling
dbh_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "dbh_normalisation.rds"))
climate_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "climate_normalisation.rds"))
ph_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "ph_normalisation.rds"))
ba_mu_sd = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), ""), "ba_normalisation.rds"))

#### Simulate growth along an environmental gradient
## Define stan model to simulate growth
gq_model = cmdstan_model("generate_growth.stan")
gq_model_classic = cmdstan_model("generate_growth_classic.stan")

## Common variables
mean_tas = mean(values(subset(climate_mean, "temperature")), na.rm = TRUE)
min_tas = min(values(subset(climate_mean, "temperature")), na.rm = TRUE)
max_tas = max(values(subset(climate_mean, "temperature")), na.rm = TRUE)

mean_pr = mean(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
min_pr = min(values(subset(climate_mean, "precipitation")), na.rm = TRUE)
max_pr = max(values(subset(climate_mean, "precipitation")), na.rm = TRUE)

mean_ph = mean(values(subset(climate_mean, "ph")), na.rm = TRUE)
min_ph = min(values(subset(climate_mean, "ph")), na.rm = TRUE)
max_ph = max(values(subset(climate_mean, "ph")), na.rm = TRUE)

n_tas = 200
n_pr = 600
n_ph = 10

n_chains = results$num_chains()

treeData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "treeData.rds"))

lower_dbh = quantile(treeData$dbh, seq(0, 1, 0.025))["2.5%"]
upper_dbh = quantile(treeData$dbh, seq(0, 1, 0.025))["97.5%"]

n_dbh_new = round(unname(round(upper_dbh - lower_dbh) + 1)/5)

## Simulate growth along a temperature gradient for many dbh
# Temperature gradient
dt_growth_tas = data.table(temperature = seq(min_tas, max_tas, length.out = n_tas), integral = numeric(n_tas))

# State-Space Model approach
stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
stanData$x_r = c((rep(mean_pr, n_tas) - stanData$pr_mu)/stanData$pr_sd, # Precip
	(dt_growth_tas[, temperature] - stanData$tas_mu)/stanData$tas_sd, # Temperature
	(rep(mean_ph, n_tas) - stanData$ph_mu)/stanData$ph_sd, # pH
	(rep(25, n_tas) - stanData$ba_mu)/stanData$ba_sd
)
stanData$n_climate_new = n_tas
stanData$n_dbh_new = n_dbh_new
stanData$n_growth = stanData$n_children # n_children is also n_growth (i.e. the number of growth intervals)!
stanData$lower_bound = unname(lower_dbh)
stanData$upper_bound = unname(upper_dbh)

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

simulatedGrowth = generate_quantities$draws("simulatedGrowth")
dim(simulatedGrowth) # n_draws x n_chains x [n_tas, n_dbh_new]

optimumTemperature = optimumPredictorValue(slope = estimated_params["tas_slope"], slope2 = estimated_params["tas_slope2"],
	scaling_mu = stanData$tas_mu, scaling_sd = stanData$tas_sd)

optimum_ind = which.min(abs(dt_growth_tas[, temperature] - optimumTemperature))

selectedData = paste0("simulatedGrowth[", optimum_ind, ",", 1:n_dbh_new, "]")
growth_dbh_fixed_temp = simulatedGrowth[, , selectedData]
growth_dbh_fixed_temp = stanData$sd_dbh*growth_dbh_fixed_temp

growth_q2.5 = apply(X = growth_dbh_fixed_temp, FUN = quantile, MARGIN = 3, probs = 0.025)
growth_q97.5 = apply(X = growth_dbh_fixed_temp, FUN = quantile, MARGIN = 3, probs = 0.975)
growth = apply(X = growth_dbh_fixed_temp, FUN = mean, MARGIN = 3)
new_dbh = seq(lower_dbh, upper_dbh, length.out = n_dbh_new)

pdf("ssm_approach.pdf", height = 10, width = 10)
plot(new_dbh, growth, xlab = "dbh", ylab = "growth", pch = 19, col = "#F4C430", ylim = c(min(growth_q2.5), max(growth_q97.5)))
polygon(c(rev(new_dbh), new_dbh), c(rev(growth_q2.5), growth_q97.5), col = "#44444422", border = NA)
dev.off()




# Classic approach
stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData_classic.rds"))
stanData$x_r = c((rep(mean_pr, n_tas) - stanData$pr_mu)/stanData$pr_sd, # Precip
	(dt_growth_tas[, temperature] - stanData$tas_mu)/stanData$tas_sd, # Temperature
	(rep(mean_ph, n_tas) - stanData$ph_mu)/stanData$ph_sd, # pH
	(rep(25, n_tas) - stanData$ba_mu)/stanData$ba_sd
)
stanData$n_climate_new = n_tas
stanData$n_dbh_new = n_dbh_new
stanData$lower_bound = unname(lower_dbh)
stanData$upper_bound = unname(upper_dbh)

generate_quantities_classic = gq_model_classic$generate_quantities(results_classic$draws(), data = stanData, parallel_chains = n_chains)

simulatedGrowth = generate_quantities_classic$draws("simulatedGrowth")
dim(simulatedGrowth) # n_draws x n_chains x [n_tas, n_dbh_new]
selectedData = paste0("simulatedGrowth[", optimum_ind, ",", 1:n_dbh_new, "]")
growth_dbh_fixed_temp = simulatedGrowth[, , selectedData]
growth_dbh_fixed_temp = stanData$sd_dbh*growth_dbh_fixed_temp

growth_q2.5_classic = apply(X = growth_dbh_fixed_temp, FUN = quantile, MARGIN = 3, probs = 0.025)
growth_q97.5_classic = apply(X = growth_dbh_fixed_temp, FUN = quantile, MARGIN = 3, probs = 0.975)
growth_classic = apply(X = growth_dbh_fixed_temp, FUN = mean, MARGIN = 3)

pdf("classic_approach.pdf", height = 10, width = 10)
plot(new_dbh, growth, xlab = "dbh", ylab = "growth", pch = 19, col = "#F4C430", ylim = c(min(growth_q2.5), max(growth_q97.5)))
polygon(c(rev(new_dbh), new_dbh), c(rev(growth_q2.5), growth_q97.5), col = "#44444422", border = NA)
dev.off()




pdf("ssm-vs-classic_approach.pdf", height = 10, width = 10)
plot(new_dbh, growth, xlab = "dbh", ylab = "growth", pch = 19, col = "#F4C430", ylim = c(min(growth_q2.5), max(growth_q97.5)))
points(new_dbh, growth_classic, pch = 19, col = "#034C4F", add = TRUE)
polygon(c(rev(new_dbh), new_dbh), c(rev(growth_q2.5), growth_q97.5), col = "#44444422", border = NA)
polygon(c(rev(new_dbh), new_dbh), c(rev(growth_q2.5_classic), growth_q97.5_classic), col = "#44444422", border = NA)
dev.off()




####! Old stuff related to the older version of generate_meanGrowth.stan (the version before January 2023).
average_growth_sim = generate_quantities$draws("average_growth_sim")
quantile_tas = data.table(q5 = numeric(n_tas), q95 = numeric(n_tas))
for (i in 1:n_tas)
	quantile_tas[i, c("q5", "q95") := as.list(sd(stanData$Yobs)*quantile(x = average_growth_sim[, , i], probs = c(0.05, 0.95)))]

# Precipitation gradient
dt_growth_pr = data.table(precipitation = seq(min_pr, max_pr, length.out = n_pr), integral = numeric(n_pr))

dt_growth_pr[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	pr = precipitation, tas = mean_tas, ph = mean_ph, basalArea = 25,
	params = estimated_params, standardised_dbh = FALSE, standardised_params = TRUE, standardised_variables = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], scaling = rbindlist(list(climate_mu_sd, ph_mu_sd, ba_mu_sd)))$value, by = precipitation]

dt_growth_pr[, integral := integral/(upper_dbh - lower_dbh)]

# Generate simulations
stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
stanData$x_r = c((dt_growth_pr[, precipitation] - stanData$pr_mu)/stanData$pr_sd, # Precip
	(rep(mean_tas, n_pr) - stanData$tas_mu)/stanData$tas_sd, # Temperature
	(rep(mean_ph, n_pr) - stanData$ph_mu)/stanData$ph_sd, # pH
	(rep(25, n_pr) - stanData$ba_mu)/stanData$ba_sd
)
stanData$n_climate_new = n_pr
stanData$lower_bound = unname(lower_dbh)
stanData$upper_bound = unname(upper_dbh)

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

average_growth_sim = generate_quantities$draws("average_growth_sim")
quantile_pr = data.table(q5 = numeric(n_pr), q95 = numeric(n_pr))
for (i in 1:n_pr)
	quantile_pr[i, c("q5", "q95") := as.list(sd(stanData$Yobs)*quantile(x = average_growth_sim[, , i], probs = c(0.05, 0.95)))]

# Basal area gradient
dt_growth_ba = data.table(basalArea = seq(min_ba, max_ba, length.out = n_ba), integral = numeric(n_ba))

dt_growth_ba[, integral := integrate(growth_fct, lower = lower_dbh, upper = upper_dbh,
	pr = mean_pr, ba = temperature, ph = mean_ph, basalArea = 25,
	params = estimated_params, standardised_dbh = FALSE, standardised_params = TRUE, standardised_variables = FALSE,
	sd_dbh = dbh_mu_sd[1, sd], scaling = rbindlist(list(climate_mu_sd, ph_mu_sd, ba_mu_sd)))$value, by = temperature]

dt_growth_ba[, integral := integral/(upper_dbh - lower_dbh)]

stanData = readRDS(paste0(tree_path, ifelse(!is.null(run), paste0(run, "_"), run), "stanData.rds"))
stanData$x_r = c((rep(mean_pr, n_ba) - stanData$pr_mu)/stanData$pr_sd, # Precip
	(dt_growth_ba[, temperature] - stanData$ba_mu)/stanData$ba_sd, # Temperature
	(rep(25, n_ba) - stanData$ba_mu)/stanData$ba_sd
)
stanData$n_climate_new = n_ba
stanData$lower_bound = unname(lower_dbh)
stanData$upper_bound = unname(upper_dbh)

generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)

average_growth_sim = generate_quantities$draws("average_growth_sim")
quantile_ba = data.table(q5 = numeric(n_ba), q95 = numeric(n_ba))
for (i in 1:n_ba)
	quantile_ba[i, c("q5", "q95") := as.list(sd(stanData$Yobs)*quantile(x = average_growth_sim[, , i], probs = c(0.05, 0.95)))]

#### Plot average growth
## In function of climate
# Common variables
y_min = min(dt_growth_tas[, min(integral)], dt_growth_pr[, integral],
	quantile_tas[, min(q5)], quantile_pr[, min(q5)])

y_max = max(dt_growth_tas[, max(integral)], dt_growth_pr[, max(integral)],
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
