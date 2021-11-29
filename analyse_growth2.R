
#### Aim of prog: Rescale growth to interprete the results

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)

#### Tool function
## Function to rescale the coefficients on the real scale (i.e., the coeffs have been computed on the normalised scale).
# Remark: Variables with _s are on the normalised scale, and variables with _r are on the real 'physical' scale
rescaleParams = function(optimal_clim_s, width_clim_niche_s, norm_clim_dt)
{
	# Get the normalising constantes from
	mu_clim = norm_clim_dt[variable == "pr", mu]
	sd_clim = norm_clim_dt[variable == "pr", sd]

	print(mu_clim)
	print(sd_clim)

	# Recompute the coeffs on real scale
	optimal_clim_r = mu_clim - sd_clim*optimal_clim_s
	width_clim_niche_r = sd_clim*width_clim_niche_s

	return (c(optimal_clim_r = optimal_clim_r, width_clim_niche_r = width_clim_niche_r))
}

real_scale_params = rescaleParams(optimal_clim_s = meanParams[["optimal_precip"]], width_clim_niche_s = meanParams[["width_precip_niche"]],
	norm_clim_dt = norm_clim_dt)

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

## Function to 
# predicted = function(dbh, precip, real_scale_params)
# {
# 	beta0 = real_scale_params["intercepts"]
# 	beta1 = real_scale_params["slopes_dbh"]
# 	beta2 = real_scale_params["slopes_precip"]
# 	beta3 = real_scale_params["quad_slopes_precip"]
# 	dbh_1 = beta0 + beta1*dbh + beta2*precip + beta3*precip^2
# 	return (dbh_1)
# }

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
real_scale_params = rescaleParams(optimal_clim_s = meanParams["optimal_precip"], width_clim_niche_s = meanParams["width_precip_niche"],
	norm_clim_dt = norm_clim_dt)


latent_1_6 = getParams(results, paste0("latentState[", 1:6, "]"))

mainFolder = "~/projects/def-dgravel/amael/postdoc/bayForDemo/BayForDemo Inventories/FR IFN/processed data/"
treeData = readRDS(paste0(mainFolder, "trees_forest_reshaped.rds"))
treeData = treeData[speciesName_sci == "Tilia_platyphyllos"]

jpeg("./latent_real.jpg", height = 1080, width = 1080, quality = 100)
plot(2000:2005, latent_1_6, pch = 19, col = "#34568B", cex = 4,
	xlab = "Year", ylab = "Diameter at breast height (scaled)")
points(x = 2000, y = (treeData[1, dbh] - norm_dbh_dt[1, mu])/norm_dbh_dt[1, sd] , pch = 19, col = "#FA7A35", cex = 2)
points(x = 2005, y = (treeData[2, dbh] - norm_dbh_dt[1, mu])/norm_dbh_dt[1, sd], pch = 19, col = "#CD212A", cex = 2)
dev.off()
