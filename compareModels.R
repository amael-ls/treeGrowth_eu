
#### Aim of prog: Compare the SSM approach to the classic approach
## Comments
# The main question is to check which model is better at predicting one annual growth. For this, we compare two likelihoods
#	of dendrochronological data D given the parameters estimated by each approach:
#		1. l_ssm(D) = [dendro | paramters_ssm, predictors]
#		2. l_classic(D) = [dendro | paramters_classic, predictors]
# For each annual growth, which one is more likely? Of course, we hope that ssm performs better, i.e., in general, we have
#	l_ssm > l_classic.


#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

# library(microbenchmark)
library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)
library(Rmpfr)
library(terra)
library(qgam)

#### Tool functions
## Source functions
source("toolFunctions.R")

## Function to extract parameters from stan object for a given chain and iter
make_theta = function(theta_ssm_all, theta_classic_all, ls_parameters, chain, iter)
{
	dim_ssm = dim(theta_ssm_all)
	dim_classic = dim(theta_classic_all)
	if (any(dim_ssm != dim_classic))
		stop("Dimensions mismatch")

	if (any((chain > dim_ssm[2]) || (chain > dim_classic[2]) || (iter >  dim_ssm[1]) || (iter > dim_classic[1])))
		stop("Out of bounds")

	theta_ssm = numeric(length(ls_parameters))
	names(theta_ssm) = ls_parameters
	theta_classic = numeric(length(ls_parameters))
	names(theta_classic) = ls_parameters

	for (current_param in ls_parameters)
	{
		theta_ssm[current_param] = theta_ssm_all[iter, chain, current_param]
		theta_classic[current_param] = theta_classic_all[iter, chain, current_param]
	}

	return (list(theta_ssm = theta_ssm, theta_classic = theta_classic))
}

## Likelihood of dendrochronological growth data given the predictors (starting dbh and environment) and known parameters theta
likelihood = function(annualGrowth_dendro, annualStarting_dbh, theta_ls, env_dt, scaling_dt, log = FALSE, output = "both")
{
	# Check-up
	# --- Dimensions
	nb_years = length(annualGrowth_dendro)
	if ((length(annualStarting_dbh) != nb_years) || (env_dt[, .N] != nb_years))
		stop("Dimensions mismatch between annualGrowth_dendro, annualStarting_dbh, and predictors")

	# --- Scalings
	if (!all(c("type", "variable", "mu", "sd") %in% names(scaling_dt)))
		stop("theta_ls must contain theta_ssm and theta_classic")
	
	setkey(scaling_dt, type, variable)

	# --- Parameters
	if (!all(c("theta_ssm", "theta_classic") %in% names(theta_ls)))
		stop("theta_ls must contain theta_ssm and theta_classic")
	
	theta_ssm = theta_ls[["theta_ssm"]]
	theta_classic = theta_ls[["theta_classic"]]
	ls_parameters = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
		"ph_slope", "ph_slope2", "competition_slope", "sigmaProc")

	if (!all(ls_parameters %in% names(theta_ssm)) || !all(ls_parameters %in% names(theta_classic)))
		stop("All the parameters must be provided in theta_ssm and theta_classic")

	# --- Predictors
	ls_predictors = c("pr", "tas", "ph", "standBasalArea")

	if (!all(ls_predictors %in% names(env_dt)))
		stop("All the predictors must be provided in env_dt")

	# --- Output option
	if (!(output %in% c("ssm", "classic", "both")))
		warning("Output option not recognised. Both (log)likelihood are returned")

	# Compute (log)likelihood
	expected_growth_meanlog_ssm = theta_ssm[["averageGrowth"]] +
		theta_ssm[["dbh_slope"]]*annualStarting_dbh/scaling_dt[.("ssm", "dbh"), sd] +
		theta_ssm[["dbh_slope2"]]*(annualStarting_dbh/scaling_dt[.("ssm", "dbh"), sd])^2 +
		theta_ssm[["pr_slope"]]*env_dt[, (pr - scaling_dt[.("ssm", "pr"), mu])/scaling_dt[.("ssm", "pr"), sd]] +
		theta_ssm[["pr_slope2"]]*env_dt[, (pr - scaling_dt[.("ssm", "pr"), mu])/scaling_dt[.("ssm", "pr"), sd]]^2 +
		theta_ssm[["tas_slope"]]*env_dt[, (tas - scaling_dt[.("ssm", "tas"), mu])/scaling_dt[.("ssm", "tas"), sd]] +
		theta_ssm[["tas_slope2"]]*env_dt[, (tas - scaling_dt[.("ssm", "tas"), mu])/scaling_dt[.("ssm", "tas"), sd]]^2 +
		theta_ssm[["ph_slope"]]*env_dt[, (ph - scaling_dt[.("ssm", "ph"), mu])/scaling_dt[.("ssm", "ph"), sd]] +
		theta_ssm[["ph_slope2"]]*env_dt[, (ph - scaling_dt[.("ssm", "ph"), mu])/scaling_dt[.("ssm", "ph"), sd]]^2 +
		theta_ssm[["competition_slope"]]*env_dt[, standBasalArea]

	if (output == "ssm")
		return (dlnorm(x = annualGrowth_dendro, meanlog = expected_growth_meanlog_ssm, sdlog = theta_ssm[["sigmaProc"]], log = log))

	expected_growth_meanlog_classic = theta_classic[["averageGrowth"]] +
		theta_classic[["dbh_slope"]]*annualStarting_dbh/scaling_dt[.("classic", "dbh"), sd] +
		theta_classic[["dbh_slope2"]]*(annualStarting_dbh/scaling_dt[.("classic", "dbh"), sd])^2 +
		theta_classic[["pr_slope"]]*env_dt[, (pr - scaling_dt[.("classic", "pr_avg"), mu])/scaling_dt[.("classic", "pr_avg"), sd]] +
		theta_classic[["pr_slope2"]]*env_dt[, (pr - scaling_dt[.("classic", "pr_avg"), mu])/scaling_dt[.("classic", "pr_avg"), sd]]^2 +
		theta_classic[["tas_slope"]]*env_dt[, (tas - scaling_dt[.("classic", "tas_avg"), mu])/scaling_dt[.("classic", "tas_avg"), sd]] +
		theta_classic[["tas_slope2"]]*env_dt[, (tas - scaling_dt[.("classic", "tas_avg"), mu])/scaling_dt[.("classic", "tas_avg"), sd]]^2 +
		theta_classic[["ph_slope"]]*env_dt[, (ph - scaling_dt[.("classic", "ph"), mu])/scaling_dt[.("classic", "ph"), sd]] +
		theta_classic[["ph_slope2"]]*env_dt[, (ph - scaling_dt[.("classic", "ph"), mu])/scaling_dt[.("classic", "ph"), sd]]^2 +
		theta_classic[["competition_slope"]]*env_dt[, standBasalArea]

	if (output == "classic")
		return (dlnorm(x = annualGrowth_dendro, meanlog = expected_growth_meanlog_classic, sdlog = theta_classic[["sigmaProc"]], log = log))

	return (list(l_ssm = dlnorm(x = annualGrowth_dendro, meanlog = expected_growth_meanlog_ssm, sdlog = theta_ssm[["sigmaProc"]], log = log),
		l_classic = dlnorm(x = annualGrowth_dendro, meanlog = expected_growth_meanlog_classic, sdlog = theta_classic[["sigmaProc"]], log = log)))
}



#? --------------------------------------------------------------------------------------------------
#* ----------------------    PART I: Compute quantiles GAM on the residuals    ----------------------
#? --------------------------------------------------------------------------------------------------
#### Load results
## Common variables
# args = c("Betula pendula", "1")			# OK, ssm better
# args = c("Fagus sylvatica", "1")			# OK, ssm better
# args = c("Picea abies", "1")				# OK, ssm better
# args = c("Pinus pinaster", "1")			# OK, ssm better
# args = c("Pinus sylvestris", "1")			# OK, ssm better
# args = c("Quercus petraea" , "1")			# OK, ssm better
args = commandArgs(trailingOnly = TRUE) # args = c("Fagus sylvatica", "1")
if (length(args) < 2)
	stop("Supply in this order the species and the run_id", call. = FALSE)

species = as.character(args[1])
run = as.integer(args[2])

tree_path = paste0("./", species, "/")
if (!dir.exists(tree_path))
	stop(paste0("Path not found for species <", species, ">."))

## Results
info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
classic = readRDS(paste0(tree_path, info_lastRun[["file"]]))

#### Load data
## Load and subset dendrochronological data
dendro = readRDS("/home/amael/project_ssm/inventories/treeRings/treeRings-climate.rds")
dendro = dendro[speciesName_sci == species]

if (dendro[, .N] == 0)
	stop("The species is not in the dendro data set")

## Add basal area (set to zero, which translates into average basal area on the real scale)
dendro[, standBasalArea := 0]

## Add a key (to make sure it is sort by plot_id, tree_id, year). Here, plot_id is redundant since it is included in plot_id
setkey(dendro, plot_id, tree_id, year)

## Load scalings
# Load scalings
dbh_scaling_ssm = readRDS(paste0(tree_path, run, "_dbh_normalisation.rds"))
ba_scaling_ssm = readRDS(paste0(tree_path, run, "_ba_normalisation.rds"))
climate_scaling_ssm = readRDS(paste0(tree_path, run, "_climate_normalisation.rds"))
ph_scaling_ssm = readRDS(paste0(tree_path, run, "_ph_normalisation.rds"))

dbh_scaling_classic = readRDS(paste0(tree_path, run, "_dbh_normalisation_classic.rds"))
ba_scaling_classic = readRDS(paste0(tree_path, run, "_ba_normalisation_classic.rds"))
climate_scaling_classic = readRDS(paste0(tree_path, run, "_climate_normalisation_classic.rds"))
ph_scaling_classic = readRDS(paste0(tree_path, run, "_ph_normalisation_classic.rds"))

scaling_ssm = rbindlist(list(dbh_scaling_ssm, ba_scaling_ssm, climate_scaling_ssm, ph_scaling_ssm))
scaling_classic = rbindlist(list(dbh_scaling_classic, ba_scaling_classic, climate_scaling_classic, ph_scaling_classic))

scaling_dt = rbindlist(list(ssm = scaling_ssm, classic = scaling_classic), idcol = "type")

## Create list of parameters
ls_parameters = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "sigmaProc")

theta_ssm_all = getParams(model_cmdstan = ssm, params_names = ls_parameters, type = "all")
theta_classic_all = getParams(model_cmdstan = classic, params_names = ls_parameters, type = "all")

#### Compute likelihood
## Common variables
n_chains = ssm$num_chains()
n_iter = ssm$metadata()$iter_sampling

ll_dt = data.table(mean_ll_ssm = numeric(n_chains * n_iter), sd_ll_ssm = numeric(n_chains * n_iter),
	mean_ll_classic = numeric(n_chains * n_iter), sd_ll_classic = numeric(n_chains * n_iter))

n_nonZero = dendro[dbh_increment_in_mm > 0, .N]

log_pred_density_ssm = matrix(data = NA, nrow = n_nonZero, ncol = n_chains*n_iter)
colnames(log_pred_density_ssm) = paste0("iter_", 1:(n_chains*n_iter))
log_pred_density_ssm = as.data.table(log_pred_density_ssm)

log_pred_density_classic = matrix(data = NA, nrow = n_nonZero, ncol = n_chains*n_iter)
colnames(log_pred_density_classic) = paste0("iter_", 1:(n_chains*n_iter))
log_pred_density_classic = as.data.table(log_pred_density_classic)

count = 1
for (chain in seq_len(n_chains))
{
	for (iter in seq_len(n_iter))
	{
		theta_ls = make_theta(theta_ssm_all, theta_classic_all, ls_parameters, chain, iter)

		like_ssm_classic = likelihood(annualGrowth_dendro = dendro[dbh_increment_in_mm > 0, dbh_increment_in_mm],
			annualStarting_dbh = dendro[dbh_increment_in_mm > 0, starting_dbh], theta_ls = theta_ls,
			env_dt = dendro[dbh_increment_in_mm > 0, .(ph, pr, tas, standBasalArea)], scaling_dt = scaling_dt, log = FALSE)

		log_pred_density_ssm[, paste0("iter_", count) := like_ssm_classic[["l_ssm"]]]
		log_pred_density_classic[, paste0("iter_", count) := like_ssm_classic[["l_classic"]]]
		
		mean_ssm_classic = lapply(X = like_ssm_classic, FUN = mean)
		sd_ssm_classic = lapply(X = like_ssm_classic, FUN = sd)
		ll_dt[count, c("mean_ll_ssm", "mean_ll_classic") := mean_ssm_classic]
		ll_dt[count, c("sd_ll_ssm", "sd_ll_classic") := sd_ssm_classic]
		count = count + 1
		if (count %% 100 == 0)
			print(paste0(round(100*count/ll_dt[, .N], 2), "% done"))
	}
}

saveRDS(ll_dt, paste0(tree_path, "mean_sd_likelihood.rds"))
saveRDS(log_pred_density_ssm, paste0(tree_path, "posterior_pred_distrib_ssm.rds"))
saveRDS(log_pred_density_classic, paste0(tree_path, "posterior_pred_distrib_classic.rds"))

if (ll_dt[, mean(mean_ll_ssm)] > ll_dt[, mean(mean_ll_classic)])
{
	print("SSM seems to be better")
} else {
	print("Classic seems to be better")
}

# The dimension of log_pred_density_ssm is n_measurements x (n_chain * n_iter), i.e., n_measurements x 3000
# In other words, each row corresponds to 3000 draws for one measured dbh increment.
# It also means that I can compute the mean for each row, i.e. in the column direction

log_pred_density_ssm = readRDS(paste0(tree_path, "posterior_pred_distrib_ssm.rds"))
log_pred_density_classic = readRDS(paste0(tree_path, "posterior_pred_distrib_classic.rds"))

start = Sys.time()
lpd_ssm = apply(X = log_pred_density_ssm, MARGIN = 2, FUN = function(x) {return(prod(mpfr(x, precBits = 10)))})
end = Sys.time()

end - start

lpd_reshape_ssm = c(lpd_ssm[[1]])
for (i in 2:length(lpd_ssm))
	lpd_reshape_ssm = rbind(lpd_reshape_ssm, lpd_ssm[[i]])

saveRDS(lpd_reshape_ssm, paste0(tree_path, "likelihood_pred_distrib_reshape_ssm.rds"))
# lpd_reshape_ssm = readRDS(paste0(tree_path, "likelihood_pred_distrib_reshape_ssm.rds"))

log_pred_density_ssm_data = log(mean(lpd_reshape_ssm))

start = Sys.time()
lpd_classic = apply(X = log_pred_density_classic, MARGIN = 2, FUN = function(x) {return(prod(mpfr(x, precBits = 10)))})
end = Sys.time()

end - start

lpd_reshape_classic = c(lpd_classic[[1]])
for (i in 2:length(lpd_classic))
	lpd_reshape_classic = rbind(lpd_reshape_classic, lpd_classic[[i]])

saveRDS(lpd_reshape_classic, paste0(tree_path, "likelihood_pred_distrib_reshape_classic.rds"))
# lpd_reshape_classic = readRDS(paste0(tree_path, "likelihood_pred_distrib_reshape_classic.rds"))

log_pred_density_classic_data = log(mean(lpd_reshape_classic))

log_pred_density_ssm_data - log_pred_density_classic_data

test_ssm = apply(X = log_pred_density_ssm, MARGIN = 1, FUN = mean)
test_classic = apply(X = log_pred_density_classic, MARGIN = 1, FUN = mean)
test = log10(test_ssm/test_classic) # This works intuitively only when both test_ssm and test_classic are positive (should be the case, proba!)


pdf(paste0(tree_path, "ll_ratio_vs_dbh.pdf"), height = 4, width = 4)
# --- Plot points
par(cex.axis = 1.1, cex.lab = 1.1, las = 1, mar = c(3.5, 4, 1, 1), mgp = c(2, 0.8, 0))
plot(dendro[dbh_increment_in_mm > 0, dbh], test, pch = 20, col = "#33223311", xlab = "Diameter (mm)",
	ylab = latex2exp::TeX(r"($ \log_{10}(ssm/classic) $)"), ylim = c(-2, max(test)))

# --- Add polygons for SSM and CLassic
polygon(x = c(0, 0, 3000, 3000), y = c(0, 100, 100, 0), col = "#E9851D33", border = FALSE)
polygon(x = c(0, 0, 3000, 3000), y = c(0, -100, -100, 0), col = "#2E77AB33", border = FALSE)

points(dendro[dbh_increment_in_mm > 0, dbh], test, pch = 20, col = "#33223355")

# --- Add zero line
abline(h = 0, lwd = 1, lty = "dashed")

# Compute qgam
data_qgam = data.table(dbh = dendro[dbh_increment_in_mm > 0, dbh], test)
fit = mqgam(test ~ s(dbh, k = 20, bs = "ad"), qu = c(0.025, 0.5, 0.975), data = data_qgam)
xSeq = data.table(dbh = seq(data_qgam[, min(dbh)], data_qgam[, max(dbh)], length.out = 5e2))
pred = qdo(obj = fit, fun = mgcv::predict.gam, newdata = xSeq, se = TRUE)

# Plot qgam
lines(xSeq[, dbh], pred[[1]]$fit, lwd = 2)
lines(xSeq[, dbh], pred[[2]]$fit, lwd = 2)
lines(xSeq[, dbh], pred[[3]]$fit, lwd = 2)

dev.off()

pdf(paste0(tree_path, "ll_ratio_vs_ph.pdf"), height = 4, width = 4)
par(cex.axis = 1.1, cex.lab = 1.1, las = 1, mar = c(3.5, 4, 1, 1), mgp = c(2, 0.8, 0))
plot(dendro[dbh_increment_in_mm > 0, ph], test, pch = 20, col = "#33223355", xlab = "pH",
	ylab = latex2exp::TeX(r"($ \log_{10}(ssm/classic) $)"), ylim = c(-2, max(test)))
abline(h = 0, lwd = 1, lty = "dashed")

data_qgam = data.table(ph = dendro[dbh_increment_in_mm > 0, ph], test)
fit = mqgam(test ~ s(ph, k = 20, bs = "ad"), qu = c(0.025, 0.5, 0.975), data = data_qgam)
xSeq = data.table(ph = seq(data_qgam[, min(ph)], data_qgam[, max(ph)], length.out = 5e2))
pred = qdo(obj = fit, fun = mgcv::predict.gam, newdata = xSeq, se = TRUE)

lines(xSeq[, ph], pred[[1]]$fit, lwd = 3)
lines(xSeq[, ph], pred[[2]]$fit, lwd = 3)
lines(xSeq[, ph], pred[[3]]$fit, lwd = 3)

dev.off()


pdf(paste0(tree_path, "ll_ratio_vs_pr.pdf"), height = 3, width = 3)
par(cex.axis = 1.1, cex.lab = 1.1, las = 1, mar = c(3.5, 4, 1, 1), mgp = c(2, 0.8, 0))
plot(dendro[dbh_increment_in_mm > 0, pr], test, pch = 20, col = "#33223311", xlab = "Precipitation (mm/year)",
	ylab = latex2exp::TeX(r"($ \log_{10}(ssm/classic) $)"), ylim = c(-2, max(test)))

# --- Add polygons for SSM and CLassic
polygon(x = c(0, 0, 3000, 3000), y = c(0, 100, 100, 0), col = "#E9851D33", border = FALSE)
polygon(x = c(0, 0, 3000, 3000), y = c(0, -100, -100, 0), col = "#2E77AB33", border = FALSE)

points(dendro[dbh_increment_in_mm > 0, pr], test, pch = 20, col = "#33223355")

# --- Add zero line
abline(h = 0, lwd = 1, lty = "dashed")

data_qgam = data.table(pr = dendro[dbh_increment_in_mm > 0, pr], test)
fit = mqgam(test ~ s(pr, k = 20, bs = "ad"), qu = c(0.025, 0.5, 0.975), data = data_qgam)
xSeq = data.table(pr = seq(data_qgam[, min(pr)], data_qgam[, max(pr)], length.out = 5e2))
pred = qdo(obj = fit, fun = mgcv::predict.gam, newdata = xSeq, se = TRUE)

lines(xSeq[, pr], pred[[1]]$fit, lwd = 2)
lines(xSeq[, pr], pred[[2]]$fit, lwd = 2)
lines(xSeq[, pr], pred[[3]]$fit, lwd = 2)

dev.off()


pdf(paste0(tree_path, "ll_ratio_vs_tas.pdf"), height = 3, width = 3)
par(cex.axis = 1.1, cex.lab = 1.1, las = 1, mar = c(3.5, 4, 1, 1), mgp = c(2, 0.8, 0))
plot(dendro[dbh_increment_in_mm > 0, tas], test, pch = 20, col = "#33223311", xlab = "Temperature (°C)",
	ylab = latex2exp::TeX(r"($ \log_{10}(ssm/classic) $)"), ylim = c(-2, max(test)))

# --- Add polygons for SSM and CLassic
polygon(x = c(-10, -10, 3000, 3000), y = c(0, 100, 100, 0), col = "#E9851D33", border = FALSE)
polygon(x = c(-10, -10, 3000, 3000), y = c(0, -100, -100, 0), col = "#2E77AB33", border = FALSE)

points(dendro[dbh_increment_in_mm > 0, tas], test, pch = 20, col = "#33223355")

# --- Add zero line
abline(h = 0, lwd = 1, lty = "dashed")

data_qgam = data.table(tas = dendro[dbh_increment_in_mm > 0, tas], test)
fit = mqgam(test ~ s(tas, k = 20, bs = "ad"), qu = c(0.025, 0.5, 0.975), data = data_qgam)
xSeq = data.table(tas = seq(data_qgam[, min(tas)], data_qgam[, max(tas)], length.out = 5e2))
pred = qdo(obj = fit, fun = mgcv::predict.gam, newdata = xSeq, se = TRUE)

lines(xSeq[, tas], pred[[1]]$fit, lwd = 2)
lines(xSeq[, tas], pred[[2]]$fit, lwd = 2)
lines(xSeq[, tas], pred[[3]]$fit, lwd = 2)

dev.off()

print(paste0(round(100*length(test[test > 0])/length(test), 2), "% of the time, SSM is better"))



#? -----------------------------------------------------------------------------------------------
#* ----------------------    PART II: Compare predictions for some trees    ----------------------
#? -----------------------------------------------------------------------------------------------
#### Explanations:
# I choose all the trees among the 6 species in common between the succesful parametrisation and the tree-ring data set
#	and then subset them such that:
#		- t0 < 2003
#		- t1 > 2003
#		- quantiles(precip_2003) < 25%
#		- quantiles(temperature_2003) > 75%
#
# With such conditions, I expect a drop around 2003 or 2004.
#	- Is this drop better detected in average by SSM approach?
#	- Is the response above consistent among dbh?

#### Load data
indices = readRDS(paste0(tree_path, "1_indices.rds"))
treeData = readRDS(paste0(tree_path, "1_treeData.rds"))
treeData = na.omit(treeData)
env = readRDS("/home/amael/project_ssm/inventories/growth/time_space.rds")$env
setkey(env, plot_id, year)

#### Subset data
## Time selection
ind = indices$indices
# ind = ind[type == "parent"][(year < 2003) & (year + nbYearsGrowth) > 2003, .(plot_id, tree_id, index_gen, index_latent_growth)]
ind = ind[(type == "parent") & (nfi %in% c("france", "germany")), .(plot_id, tree_id, index_gen, index_latent_growth)]
# ind = ind[, .(plot_id, tree_id, index_gen, index_latent_growth, nbIntervalGrowth)]


selectedTrees = treeData[.(ind[, .(plot_id, tree_id)])]
ls_years = selectedTrees[, min(year)]:selectedTrees[, max(year + deltaYear)]
coords = unique(selectedTrees[, .(plot_id, x, y)])

if (length(unique(coords[, plot_id])) != coords[, .N])
	stop("I assumed hereafter that plot_id is equivalent to coordinates, i.e., that they are unique!")

env = env[.(coords[, plot_id], 2003), .(plot_id, pr, tas)]

## Compute quantiles
pr_quantiles = climQuantiles(clim_var = "pr", years = ls_years, coords = coords)
tas_quantiles = climQuantiles(clim_var = "tas", years = ls_years, coords = coords)

## Merge env and the quantiles
env = merge.data.table(x = env, y = pr_quantiles[, .(plot_id, q05)], by = "plot_id")
setnames(env, old = "q05", new = "pr_05")

env = merge.data.table(x = env, y = tas_quantiles[, .(plot_id, q95)], by = "plot_id")
setnames(env, old = "q95", new = "tas_95")

## Quantile selection
selected_plots = env[(tas > tas_95) & (pr < pr_05), plot_id]

## Final subset
treeData = treeData[.(selected_plots)]
selectedTrees = unique(treeData[, .(plot_id, tree_id)])
setkey(ind, plot_id, tree_id)
ind = ind[.(selectedTrees)]

treeData = merge.data.table(x = treeData, y = ind, by = c("plot_id", "tree_id"))
treeData[, index_gen_end := index_gen + deltaYear - 1]
treeData[, index_latent_growth_end := index_latent_growth + deltaYear - 1]

#### Make figure
dbh_ls = growth_timeSeries(species = species, run = run, selected_plot_id = treeData[, plot_id], init_dbh = treeData[, dbh],
	nbYearsGrowth_new = treeData[, deltaYear], year_start_indiv = treeData[, year], caption = TRUE, extension = c("pdf", "tex"))

saveRDS(dbh_ls, paste0(tree_path, "dbh_ls.rds"))



####! CRASH TEST ZONE ------------------------------------------------------------------------
#### Extract draws
## SSM
ssm = readRDS("./Fagus sylvatica/growth-run=1-2023-01-12_02h02_de-fr-sw_12000_main.rds")
aa = ssm$draws("latent_growth")
dim(aa)

results_dt = data.table(sd_timeSeries = numeric(treeData[, .N]), sd_timeSeries_ss = numeric(treeData[, .N]))

for (i in seq_along(results))
{
	ss = dbh_scaling_ssm[variable == "dbh", sd]*aa[, , treeData[i, index_latent_growth]:treeData[i, index_latent_growth_end]]
	meanG = apply(X = ss, FUN = mean, MARGIN = 3)
	results_dt[i, sd_timeSeries := sd(meanG)]
	results_dt[i, sd_timeSeries_ss := sd(ss)]

	if (i %% 100 == 0)
		print(paste0(round(100*i/treeData[, .N]), "% done"))
}

pdf("Histogram_variance.pdf", height = 3, width = 3)
par(cex.axis = 1.25, cex.lab = 1.25, las = 1, mar = c(3.5, 4, 1, 1), mgp = c(2, 0.8, 0))
hist(results_dt[, log(sd_timeSeries)], xlab = "log(std. dev.)", ylab = "", main = "")
dev.off()


hist(results_dt[, sd_timeSeries], xlab = "log(std. dev.)", ylab = "", main = "", prob = TRUE)
curve(dlnorm(x, meanlog = mean(results_dt[, log(sd_timeSeries)]), sdlog = sd(results_dt[, log(sd_timeSeries)])), lwd = 2, add = TRUE)
meanlog = mean(results_dt[, log(sd_timeSeries)])
sdlog = sd(results_dt[, log(sd_timeSeries)])

exp(meanlog + sdlog^2/2)

pdf(paste0(tree_path, "time_series_ssm_1.pdf"), height = 2.472136, width = 4) # 2/3 golden ratio
plot(treeData[i, year]:treeData[i, year + deltaYear - 1], meanG, type = "l", xlab = "Year",
	ylab = "Predicted growth", las = 1)
points(treeData[i, year]:treeData[i, year + deltaYear - 1], meanG, pch = 19)
dev.off()



# * ----------------------------------------------------------------



# simulated_dbh_ssm = dbh_ls[["simulated_dbh_ssm"]]
# simulated_dbh_classic = dbh_ls[["simulated_dbh_classic"]]
# simulated_dbh_classic_avg = dbh_ls[["simulated_dbh_classic_avg"]]

# years = dbh_ls[["years"]]

# simulated_dbh_ssm_avg = apply(X = simulated_dbh_ssm, FUN = mean, MARGIN = 3)
# simulated_dbh_classic_avg = apply(X = simulated_dbh_classic, FUN = mean, MARGIN = 3)

# lm_ssm = lm(simulated_dbh_ssm_avg ~ years)
# lm_classic = lm(simulated_dbh_classic_avg ~ years)

# plot(years, simulated_dbh_ssm_avg, type = "l")
# points(years, simulated_dbh_ssm_avg, pch = 19)
# abline(lm_ssm, lwd = 2, col = "#CD211A")

# plot(years, simulated_dbh_classic_avg, , type = "l")
# points(years, simulated_dbh_classic_avg, pch = 19)
# abline(lm_classic, lwd = 2, col = "#CD211A")


# reload = readRDS(paste0(tree_path, "1_treeData.rds"))
# reload = na.omit(reload)
# ii = reload[.(treeData[2, .(plot_id, tree_id)]), which = TRUE]
# reload[ii]
# dbh_reload = reload[.(treeData[2, .(plot_id, tree_id)]), dbh]

# bb = ssm$draws("latent_dbh_parents")
# ss_dbh = dbh_scaling_ssm[variable == "dbh", sd]*bb[, , ii[1]]
# lazyPosterior(ss_dbh)



# * ----------------------------------------------------------------



# simulated_dbh_ssm = generate_quantities_ssm$draws("current_dbh")

# simulatedGrowth_avg_ssm = dbh_scaling_ssm[variable == "dbh", sd]*generate_quantities_ssm$draws("simulatedGrowth_avg")
# simulatedGrowth_avg_classic = dbh_scaling_ssm[variable == "dbh", sd]*generate_quantities_classic$draws("simulatedGrowth_avg")
# simulatedGrowth_avg_clim_avg = dbh_scaling_ssm[variable == "dbh", sd]*generate_quantities_classic$draws("simulatedGrowth_avg_clim_avg")

# simulatedGrowth_ssm = generate_quantities_ssm$draws("simulatedGrowth")

# aa = simulated_dbh_ssm[, , 1] + simulatedGrowth_ssm[, , 1] - simulated_dbh_ssm[, , 2]


# simulatedGrowth_ssm = dbh_scaling_ssm[variable == "dbh", sd]*simulatedGrowth_ssm
# simulatedGrowth_ssm_qt = apply(X = simulatedGrowth_ssm, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.5, 0.95))

# simulatedGrowth_avg_ssm_qt = apply(X = simulatedGrowth_avg_ssm, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.5, 0.95))

# simulatedGrowth_ssm_mean = apply(X = simulatedGrowth_ssm, MARGIN = 3, FUN = mean)

# simulatedGrowth_avg_ssm_mean = apply(X = simulatedGrowth_avg_ssm, MARGIN = 3, FUN = mean)

# lims = c(min(simulatedGrowth_ssm_qt["5%", ]), max(simulatedGrowth_ssm_qt["95%", ]))
# lims = c(2.2, 2.8)

# plot(year_start:year_end, simulatedGrowth_ssm_qt["50%",], type = "l", xlab = "Year", ylab = "Growth", lwd = 3, col = "#F4C430",
# 	ylim = lims)
# polygon(c(rev(year_start:year_end), year_start:year_end), c(rev(simulatedGrowth_ssm_qt["5%",]), simulatedGrowth_ssm_qt["95%",]),
# 	col = "#F4C43022", border = NA)

# lines(year_start:year_end, simulatedGrowth_avg_ssm_qt["50%",], type = "l", xlab = "Year", ylab = "Growth", lwd = 3, col = "#3060f4")
# polygon(c(rev(year_start:year_end), year_start:year_end), c(rev(simulatedGrowth_avg_ssm_qt["5%",]), simulatedGrowth_avg_ssm_qt["95%",]),
# 	col = "#3060f422", border = NA)

# lines(year_start:year_end, simulatedGrowth_ssm_mean, col = "#4F0603", lwd = 3, lty = 2)

# lines(year_start:year_end, simulatedGrowth_avg_ssm_mean, col = "#BB7799", lwd = 3, lty = 3)

####! END CRASH TEST ZONE ------------------------------------------------------------------------
