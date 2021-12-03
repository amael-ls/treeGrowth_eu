
#### Aim of prog: Rescale growth to interprete the results

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)
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
	return (dt[.N, file])
}

#### Load results
results = readRDS("Tilia_platyphyllos/growth-2021-11-29_05h17.rds")

isPrecip_normalised = TRUE
isDBH_normalised = FALSE

if (isPrecip_normalised)
{
	print("Precip parameter must be transformed when working on the real precip scale")
	norm_clim_dt = readRDS("Tilia_platyphyllos/climate_normalisation.rds")
}

if (isDBH_normalised)
{
	print("DBH must be transformed when working on the real DBH scale")
	norm_dbh_dt = readRDS("Tilia_platyphyllos/dbh_normalisation.rds")
}

#### Rescale parameters
paramsNames = c("potentialMaxGrowth", "power_dbh", "optimal_precip", "width_precip_niche", "processError")
meanParams = getParams(results, paramsNames)
real_scale_params = rescaleParams(optimal_clim_s = meanParams["optimal_precip"],
	width_clim_niche_s = meanParams["width_precip_niche"], norm_clim_dt = norm_clim_dt)

latent_1_6 = getParams(results, paste0("latent_dbh[", 1:6, "]"))
x = 2000:2005
test = lm(latent_1_6 ~ x)
summary(test)

mainFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
treeData = treeData[speciesName_sci == "Tilia_platyphyllos"]

jpeg("./latent_real_2ndOption.jpg", height = 1080, width = 1080, quality = 100)
plot(2000:2005, latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)")
points(x = 2000, y = treeData[1, dbh], pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = treeData[2, dbh], pch = 19, col = "#CD212A", cex = 2)
abline(test)
dev.off()

pdf("./growth_dbh.pdf", height = 7, width = 7)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.5, 0), tck = -0.015)
curve(growth(x, 0, params = meanParams, norm_clim_dt, isScaled = TRUE), from = 50, to = 1000,
	lwd = 2, col = "#34568B", xlab = "dbh", ylab = "Growth (in mm/yr)")
dev.off()

####! CRASH TEST ZONE
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