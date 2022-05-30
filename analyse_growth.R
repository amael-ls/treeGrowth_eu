
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

#### Tool functions
source("./toolFunctions.R")

#### Common variables
## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

n_runs = 4 # Number of runs used in growth_subsample.R
threshold_indiv = 8000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2022/05/17") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(speciesName_sci, multiRun, threshold_time, lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

error_ls = vector(mode = "list", length = infoSpecies[, .N])
names(error_ls) = infoSpecies[, speciesName_sci]

correl_ls = vector(mode = "list", length = infoSpecies[, .N])
names(correl_ls) = infoSpecies[, speciesName_sci]

posterior_ls = vector(mode = "list", length = infoSpecies[, .N])
names(posterior_ls) = infoSpecies[, speciesName_sci]

params_dt = data.table(parameters = c("averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2",
	"tas_slope", "tas_slope2", "ph_slope", "ph_slope2", "competition_slope"),
	priors = c(dnorm, dgamma, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm),
	arg1 = c(-4, 1/10, 0, 0, -1, 0, -1, 0, -1, 0),
	arg2 = c(2, 1/10, 5, 5, 0.75, 5, 0.75, 5, 0.75, 5),
	title = c("Average growth (mean)", "Random effect (sd)", "Dbh slope", "Precipitation slope", "Precipitation slope (quadratic term)",
		"Temperature slope", "Temperature slope (quadratic term)", "Soil acidity slope (pH)", "Soil acidity slope (pH, quadratic term)",
		"Competition slope"),
	expand_bounds = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

setkey(params_dt, parameters)

#### For loop on processed species, to plot posteriors of the main parameters (errors, intercept, and slopes)
for (species in infoSpecies[, speciesName_sci])
{
	multi = infoSpecies[species, multiRun]
	ls_nfi = unlist(stri_split(infoSpecies[species, ls_nfi], regex = ", "))
	summary_dt = centralised_fct(species, multi, n_runs, ls_nfi, params_dt, run = if (multi) NULL else 1)
	error_ls[[species]] = summary_dt[["error_dt"]]
	correl_ls[[species]] = summary_dt[["correl_energy"]]
	posterior_ls[[species]] = summary_dt[["posteriorSim"]]
}

error_dt = rbindlist(error_ls, idcol = "speciesName_sci")
correl_dt = rbindlist(correl_ls, idcol = "speciesName_sci")

saveRDS(error_dt, "./error_species.rds")
saveRDS(correl_dt, "./correlation_energy_species.rds")

plot_correl_error(error_dt, correl_dt, threshold_correl = 0.2, rm_correl = "lp__")

selected_indiv = data.table(
	speciesName_sci = infoSpecies[, speciesName_sci],
	plot_id = c("france_738952", "france_804256", "france_873749"),
	tree_id = c(1, 3, 5)
)

setkey(selected_indiv, speciesName_sci)

# computeGrowth(dt = treeData, byCols = c("plot_id", "tree_id"))
abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries.pdf",
	plot_id = selected_indiv["Abies grandis", plot_id], tree_id = selected_indiv["Abies grandis", tree_id],
	highlight_threshold_growth = 45)

abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries.pdf",
	plot_id = "france_725338", tree_id = 16,
	highlight_threshold_growth = 50)

abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries2.pdf",
	plot_id = "france_818095", tree_id = 4,
	highlight_threshold_growth = 50)

acerOpalus = dbh_timeSeries(posterior_ls[["Acer opalus"]], plotMean = TRUE, filename = "timeSeries.pdf", data = treeData,
	plot_id = selected_indiv["Acer opalus", plot_id], tree_id = selected_indiv["Acer opalus", tree_id],
	highlight_threshold_growth = 45)

tiliaPlat = dbh_timeSeries(posterior_ls[["Tilia platyphyllos"]], plotMean = TRUE, filename = "timeSeries.pdf", data = treeData,
	plot_id = selected_indiv["Tilia platyphyllos", plot_id], tree_id = selected_indiv["Tilia platyphyllos", tree_id],
	highlight_threshold_growth = 45)

#* --------------------------------------------------------------------------------------------
#! ********************************************************************************************
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================    WHAT FOLLOWS ARE OLDER STUFF, SOME ARE GOOD    ================ ODOT#
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! ********************************************************************************************
#* --------------------------------------------------------------------------------------------

#! TO REMOVE AFTER !!!!!!!!!!!!!!!
results = readRDS("Abies grandis/growth-run=1-2022-05-18_00h19.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
stanData = readRDS("Abies grandis/1_stanData.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
indices = readRDS("Abies grandis/1_indices.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
treeData = readRDS("Abies grandis/1_treeData.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
divergences = which(results$sampler_diagnostics()[, , "divergent__"] == 1) #! TO REMOVE AFTER !!!!!!!!!!!!!!!
sd_dbh = readRDS("Abies grandis/1_dbh_normalisation.rds")[, sd] #! TO REMOVE AFTER !!!!!!!!!!!!!!!
clim_scaling = readRDS("Abies grandis/1_climate_normalisation.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
ph_scaling = readRDS("Abies grandis/1_ph_normalisation.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
ba_scaling = readRDS("Abies grandis/1_ba_normalisation.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!

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
dim(generate_quantities$draws()) # iter_sampling * n_chains * (2*n_obs + n_latentGrowth + n_hiddenState + 2), +2 from simplex

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
# abline(v = stanData$n_indiv, lwd = 3, col = "#CD212A")
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

####! CRASH TEST ZONE
yearly_latent_dbh_array = generate_quantities$draws("yearly_latent_dbh") # iter_sampling * n_chains * n_hiddenState
yearly_latent_dbh_array = sd_dbh*yearly_latent_dbh_array

## With first tree
temporary_array = yearly_latent_dbh_array[, , indices[1, index_gen]:indices[2, index_gen]]
dimnames(temporary_array)$variable = indices[1, year]:indices[2, year]
np = bayesplot::nuts_params(results)
str(np)
levels(np$Parameter)

bayesplot::color_scheme_set("darkgray")
div_style = bayesplot::parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
pdf("test.pdf", height = 11.25, width = 20)
bayesplot::mcmc_parcoord(temporary_array, size = 0.25, alpha = 0.1, np = np, np_style = div_style)
dev.off()

## With tree parent number 812, corresponds to indices[1729]
temporary_array = yearly_latent_dbh_array[, , indices[1729, index_gen]:indices[1730, index_gen]]
dimnames(temporary_array)$variable = indices[1729, year]:indices[1730, year]
np = bayesplot::nuts_params(results)
str(np)
levels(np$Parameter)

bayesplot::color_scheme_set("darkgray")
div_style = bayesplot::parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
pdf("test_1729.pdf", height = 11.25, width = 20)
bayesplot::mcmc_parcoord(temporary_array, size = 0.25, alpha = 0.1, np = np, np_style = div_style)
dev.off()

draws = temporary_array
data = treeData[1:2]
data = treeData[1729:1730]

# r$> selectedIndiv
# [1] 1

# r$> selectedPlot
# [1] "france_738952"

ll = dbh_timeSeries(posteriorSim, plotMean = TRUE, filename = "test.pdf", plot_id = "france_738952", tree_id = 1,
	highlight_threshold_growth = 45)



avg_G = results$draws("averageGrowth")
plotEffect = results$draws("plotEffect")

sd_dbh*exp(quantile(avg_G))
quantile(plotEffect)

params_dt = data.table(ls_params = c("averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope",
	"tas_slope2", "ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), value = 0)

params_dt[, isPredictor := FALSE]
params_dt[4:10, isPredictor := TRUE]

setkey(params_dt, ls_params)

for (param in params_dt[, ls_params]) #! A little stupid, there is the getParams function...
	params_dt[param, value := mean(results$draws(param))]


setkey(clim_scaling, variable)
setkey(ph_scaling, variable)
setkey(ba_scaling, variable)

params_dt[stri_detect(ls_params, regex = "pr_slo"), mu := clim_scaling["pr", mu]]
params_dt[stri_detect(ls_params, regex = "pr_slo"), sd := clim_scaling["pr", sd]]

params_dt[stri_detect(ls_params, regex = "tas_slo"), mu := clim_scaling["tas", mu]]
params_dt[stri_detect(ls_params, regex = "tas_slo"), sd := clim_scaling["tas", sd]]

params_dt[stri_detect(ls_params, regex = "ph_slo"), mu := ph_scaling["ph", mu]]
params_dt[stri_detect(ls_params, regex = "ph_slo"), sd := ph_scaling["ph", sd]]

params_dt[stri_detect(ls_params, regex = "competit"), mu := ba_scaling["standBasalArea_interp", mu]]
params_dt[stri_detect(ls_params, regex = "competit"), sd := ba_scaling["standBasalArea_interp", sd]]

# rescaleParams(intercept = params_dt["averageGrowth_mu", value], slope_dbh = params_dt["dbh_slope", value],
# 	slope_predictors = params_dt[(isPredictor), value], sd_dbh = sd_dbh,
# 	mu_predictors = params_dt[(isPredictor), mu], sd_predictors = params_dt[(isPredictor), sd])

scaling = rbindlist(list(clim_scaling, ph_scaling, ba_scaling))
setkey(scaling, variable)

params = rescaleParams(params = getParams(results, c("averageGrowth_mu", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope",
		"tas_slope2", "ph_slope", "ph_slope2", "competition_slope")),
	sd_dbh = sd_dbh,
	mu_predictors = c(pr = scaling["pr", mu], tas = scaling["tas", mu], ph = scaling["ph", mu],
		basalArea = scaling["standBasalArea_interp", mu]),
	sd_predictors = c(pr = scaling["pr", sd], tas = scaling["tas", sd], ph = scaling["ph", sd],
		basalArea = scaling["standBasalArea_interp", sd]))

growth_fct(250, 750, 12, 5, 25, params, sd_dbh, rescaled = TRUE)
growth_fct(500, 500, 12, 5, 25, params, sd_dbh, rescaled = TRUE)

