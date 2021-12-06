
#### Aim of prog: Rescale growth to interprete the results

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(bayesplot)
library(cmdstanr)
library(stringi)
library(DHARMa)
library(terra)

#### Tool function
## Function to rescale the coefficients on the real scale (i.e., the coeffs have been computed on the normalised scale).
# Remark: Variables with _s are on the normalised scale, and variables with _r are on the real 'physical' scale
rescaleParams = function(optimal_clim_s, width_clim_niche_s, norm_clim_dt)
{
	# Get the normalising constantes from
	mu_clim = norm_clim_dt[variable == "pr", mu]
	sd_clim = norm_clim_dt[variable == "pr", sd]

	# Recompute the coeffs on real scale
	optimal_clim_r = mu_clim - sd_clim*optimal_clim_s
	width_clim_niche_r = sd_clim*width_clim_niche_s

	return (c(optimal_clim_r = optimal_clim_r, width_clim_niche_r = width_clim_niche_r))
}

## Get fixed values parameters
getParams = function(model_cmdstan, params_names, type = "mean")
{
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

## Growth function (expected growth)
growth = function(dbh, precip, params, norm_clim_dt, isScaled = FALSE)
{
	# Get the normalising constantes
	mu_clim = norm_clim_dt[variable == "pr", mu]
	sd_clim = norm_clim_dt[variable == "pr", sd]

	# Scale the climate
	if (isScaled)
		clim = precip
	if(!isScaled)
		clim = (precip - mu_clim)/sd_clim

	potentialMaxGrowth = params[["potentialMaxGrowth"]]
	power_dbh = params[["power_dbh"]]
	optimal_clim = params[["optimal_precip"]]
	width_clim_niche = params[["width_precip_niche"]]
	G = potentialMaxGrowth * dbh^power_dbh *
		exp(-(clim - optimal_clim)^2/width_clim_niche^2)
	return (G)
}

## Dbh function (expected dbh at time t + 1 given climate and dbh at time t)
nextDBH = function(currentDBH, precip, params, norm_clim_dt, isScaled = FALSE)
	return (currentDBH + growth(currentDBH, precip, params, norm_clim_dt, isScaled))

## Get name of the last run
getLastRun = function(path, begin = "growth-", extension = ".rds", format = "ymd")
{
	if (format != "ymd")
		stop("Format is not recognised. For now only year-month-day alias ymd works")
	
	ls_files = list.files(path = path, pattern = paste0("^", begin, ".*", extension, "$"))
	ls_files_split = stri_split(stri_sub(ls_files, from = stri_locate(ls_files, regex = begin)[, "end"] + 1),
			regex = "-", simplify = TRUE)
	n = length(ls_files)
	
	if (format == "ymd") # year month day
		dt = data.table(file = ls_files, year = ls_files_split[, 1], month = ls_files_split[, 2], day = ls_files_split[, 3])

	dt[stri_detect(str = day, regex = extension), day := stri_sub(str = day, to = stri_locate_first(str = day, regex = "_")[,"start"] - 1)]
	setorder(dt, year, month, day)
	return (list(file = dt[.N, file], time_ended = paste(dt[.N, year], dt[.N, month], dt[.N, day], sep = "-")))
}

#### Read data
## Common variables
species = "Tilia_platyphyllos"
path = paste0("./", species, "/")
info_lastRun = getLastRun(path)
lastRun = info_lastRun[["file"]]
time_ended = info_lastRun[["time_ended"]]

isPrecip_normalised = TRUE
isDBH_normalised = TRUE

## Load results
results = readRDS(paste0(path, lastRun))

if (isPrecip_normalised)
{
	print("Precip parameter must be transformed when working on the real precip scale")
	norm_clim_dt = readRDS(paste0(path, "climate_normalisation.rds"))
}

if (isDBH_normalised)
{
	print("DBH must be transformed when working on the real DBH scale")
	norm_dbh_dt = readRDS(paste0(path, "dbh_normalisation.rds"))
}

#### Get parameters
paramsNames = c("potentialMaxGrowth", "power_dbh", "optimal_precip", "width_precip_niche", "processError")
meanParams = getParams(results, paramsNames)
processError = meanParams[["processError"]]

#### Residuals
## Load data
treeFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
treeData = readRDS(paste0(treeFolder, "trees_forest_reshaped.rds"))
treeData = treeData[speciesName_sci == species]

## Get dbh and time
dbh_start = treeData[treeData[, .I[which.min(year)], by = .(tree_id, pointInventory_id)][, V1], dbh]
dbh_end = treeData[treeData[, .I[which.max(year)], by = .(tree_id, pointInventory_id)][, V1], dbh]

t_start = treeData[treeData[, .I[which.min(year)], by = .(tree_id, pointInventory_id)][, V1], year]
t_end = treeData[treeData[, .I[which.max(year)], by = .(tree_id, pointInventory_id)][, V1], year]
delta_t = t_end - t_start

## Climate
climFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/climateData/Chelsa/yearlyAverage/"
climate = readRDS(paste0(climFolder, "FR_reshaped_climate.rds"))

## indices
indices = readRDS(paste0(treeFolder, species, "_indices.rds"))
indices = unique(indices[, .(tree_id, pointInventory_id, index_clim_start, index_clim_end)])

## Create simulations
n_rep = 500
dt = data.table(rep_dbh_end = rep(dbh_end, each = n_rep), sampled = numeric(n_rep * length(dbh_end)))
if (isDBH_normalised)
	dt[, rep_dbh_end := rep_dbh_end/norm_dbh_dt[, sd]]

for (i in 1:length(dbh_end))
{
	currentDbh = dbh_start[i]/ifelse(isDBH_normalised, norm_dbh_dt[, sd], 1)
	currentClim_index = indices[i, index_clim_start]
	for (yr in 1:delta_t[i])
	{
		precip = climate[currentClim_index, pr]
		currentDbh = rnorm(n_rep, rnorm(n_rep, nextDBH(currentDbh, precip, meanParams, norm_clim_dt), processError), 0.006)
		currentClim_index = currentClim_index + 1
	}
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := currentDbh]
	if (i %% 100 == 0)
		print(paste0(round(i*100/length(dbh_end), 2), "% done"))
}

saveRDS(dt, "testForResiduals.rds")

## Reshape simulations into a matrix length(dbh) x n_rep
sims = matrix(data = dt[, sampled], nrow = n_rep, ncol = length(dbh_end)) # each column is for one data point (the way I organised the dt)
sims = t(sims) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = dbh_end/ifelse(isDBH_normalised, norm_dbh_dt[, sd], 1))
#, fittedPredictedResponse = avg(dbh, alpha0_ptw_med, alpha1_ptw_med))

plot(forDharma)
dev.off()

#### Plot posterior distributions and traces of parameters
# Folder to save plots
figurePath = paste0(path, time_ended, "/")
if (!dir.exists(figurePath))
	dir.create(figurePath)

# Intercept
plot_title = ggplot2::ggtitle("Posterior distribution potentialMaxGrowth", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("potentialMaxGrowth"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "pmg-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for potentialMaxGrowth")
temp_plot = mcmc_trace(results$draws("potentialMaxGrowth")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "pmg-traces.pdf"), plot = temp_plot, device = "pdf")

# Slopes_dbh
plot_title = ggplot2::ggtitle("Posterior distribution power_dbh", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("power_dbh"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "power_dbh-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for power_dbh")
temp_plot = mcmc_trace(results$draws("power_dbh")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "power_dbh-traces.pdf"), plot = temp_plot, device = "pdf")

# Slopes_precip
plot_title = ggplot2::ggtitle("Posterior distribution optimal_precip", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("optimal_precip"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "optimal_precip-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for optimal_precip")
temp_plot = mcmc_trace(results$draws("optimal_precip")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "optimal_precip-traces.pdf"), plot = temp_plot, device = "pdf")

# Quad_slopes_precip
plot_title = ggplot2::ggtitle("Posterior distribution width_precip_niche", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("width_precip_niche"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "width_precip_niche-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for width_precip_niche")
temp_plot = mcmc_trace(results$draws("width_precip_niche")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "width_precip_niche-traces.pdf"), plot = temp_plot, device = "pdf")

# processError
plot_title = ggplot2::ggtitle("Posterior distribution processError", "with medians and 80% intervals")
temp_plot = mcmc_areas(results$draws("processError"), prob = 0.8) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "processError-distrib.pdf"), plot = temp_plot, device = "pdf")

plot_title = ggplot2::ggtitle("Traces for processError")
temp_plot = mcmc_trace(results$draws("processError")) + plot_title
ggplot2::ggsave(filename = paste0(figurePath, "processError-traces.pdf"), plot = temp_plot, device = "pdf")

# # measureError
# plot_title = ggplot2::ggtitle("Posterior distribution measureError", "with medians and 80% intervals")
# temp_plot = mcmc_areas(results$draws("measureError"), prob = 0.8) + plot_title
# ggplot2::ggsave(filename = paste0(figurePath, "measureError-distrib.pdf"), plot = temp_plot, device = "pdf")

# plot_title = ggplot2::ggtitle("Traces for measureError")
# temp_plot = mcmc_trace(results$draws("measureError")) + plot_title
# ggplot2::ggsave(filename = paste0(figurePath, "measureError-traces.pdf"), plot = temp_plot, device = "pdf")

for (i in c(1:6, 12:14))
{
	plot_title = ggplot2::ggtitle(paste0("Posterior distribution latent_dbh[", i, "]"),
		"with medians and 80% intervals")
	temp_plot = mcmc_areas(results$draws(paste0("latent_dbh[", i, "]")), prob = 0.8) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latent_dbh_", i, "-distrib.pdf"), plot = temp_plot, device = "pdf")

	plot_title = ggplot2::ggtitle(paste0("Traces for latent_dbh[", i, "]"))
	temp_plot = mcmc_trace(results$draws(paste0("latent_dbh[", i, "]"), inc_warmup = TRUE)) + plot_title
	ggplot2::ggsave(filename = paste0(figurePath, "latent_dbh_", i, "-traces.pdf"), plot = temp_plot, device = "pdf")
}

####! CRASH TEST ZONE
latent_1_6 = getParams(results, paste0("latent_dbh[", 1:6, "]"))
x = 2000:2005

pdf("./latent_real_2ndOption.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(2000:2005, norm_dbh_dt[, sd]*latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)")
points(x = 2000, y = treeData[1, dbh], pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = treeData[2, dbh], pch = 19, col = "#CD212A", cex = 2)
dev.off()

pdf("./growth_dbh.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
curve(growth(x, 0, params = meanParams, norm_clim_dt, isScaled = TRUE), from = 50, to = 1000,
	lwd = 2, col = "#34568B", xlab = "dbh", ylab = "Growth (in mm/yr)")
dev.off()

pdf("./growth_precip.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
curve(growth(250/norm_dbh_dt[, sd], x, params = meanParams, norm_clim_dt, isScaled = FALSE), from = 650, to = 1500,
	lwd = 2, col = "#34568B", xlab = "Precipitation (mm/yr)", ylab = "Growth (in mm/yr)")
dev.off()

# What follows works only for the French data that have a particular structure. For the general case, use indices
dbh0 = treeData[seq(1, .N - 1, by = 2), dbh]
dbh1 = treeData[seq(2, .N, by = 2), dbh]
g_5_years = (dbh1 - dbh0)/5

# jpeg("test.jpg", width = 1080, height = 1080, quality = 100)
pdf("test.pdf", width = 7, height = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
plot(dbh0, g_5_years, pch = 19, col = "#FA7A354A")
dev.off()

####! END CRASH TEST ZONE