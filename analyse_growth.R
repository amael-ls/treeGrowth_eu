
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
threshold_time = as.Date("2022/07/25") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(speciesName_sci, multiRun, threshold_time, lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

error_ls = vector(mode = "list", length = infoSpecies[, .N])
names(error_ls) = infoSpecies[, speciesName_sci]

correl_ls = vector(mode = "list", length = infoSpecies[, .N])
names(correl_ls) = infoSpecies[, speciesName_sci]

posterior_ls = vector(mode = "list", length = infoSpecies[, .N])
names(posterior_ls) = infoSpecies[, speciesName_sci]

ls_files = character(length = infoSpecies[, .N])
names(ls_files) = infoSpecies[, speciesName_sci]

params_dt = data.table(parameters = c("averageGrowth", "dbh_slope", "pr_slope", "pr_slope2",
	"tas_slope", "tas_slope2", "ph_slope", "ph_slope2", "competition_slope"),
	priors = c(dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm),
	arg1 = c(-4, 0, 0, 0, 0, 0, 0, 0, 0),
	arg2 = c(10, 5, 5, 5, 5, 5, 5, 5, 5),
	title = c("Average growth (mean)", "Dbh slope", "Precipitation slope", "Precipitation slope (quadratic term)",
		"Temperature slope", "Temperature slope (quadratic term)", "Soil acidity slope (pH)", "Soil acidity slope (pH, quadratic term)",
		"Competition slope"),
	expand_bounds = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

setkey(params_dt, parameters)

#### For loop on processed species, to plot posteriors of the main parameters (errors, intercept, and slopes)
for (species in infoSpecies[, speciesName_sci])
{
	multi = infoSpecies[species, multiRun]
	ls_nfi = unlist(stri_split(infoSpecies[species, ls_nfi], regex = ", "))
	summary_dt = centralised_fct(species, multi, n_runs, ls_nfi, params_dt, run = if (multi) NULL else 1, simulatePosterior = TRUE)
	error_ls[[species]] = summary_dt[["error_dt"]]
	correl_ls[[species]] = summary_dt[["correl_energy"]]
	posterior_ls[[species]] = summary_dt[["posteriorSim"]]
	ls_files[[species]] = summary_dt[["fileResults"]]
}

error_dt = rbindlist(error_ls, idcol = "speciesName_sci")
correl_dt = rbindlist(correl_ls, idcol = "speciesName_sci")

saveRDS(error_dt, "./error_species.rds")
saveRDS(correl_dt, "./correlation_energy_species.rds")
saveRDS(posterior_ls, "./posterior_ls.rds")

plot_correl_error(error_dt, correl_dt, threshold_correl = 0.2, rm_correl = "lp__")

# selected_indiv = data.table(speciesName_sci = infoSpecies[, speciesName_sci],
# 	plot_id = c("france_738952", "france_804256", "france_873749"), tree_id = c(1, 3, 5))

# setkey(selected_indiv, speciesName_sci)

# computeGrowth(dt = treeData, byCols = c("plot_id", "tree_id"))
# abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries.pdf",
# 	plot_id = selected_indiv["Abies grandis", plot_id], tree_id = selected_indiv["Abies grandis", tree_id],
# 	highlight_threshold_growth = 45)

#### Few trajectory plots
## Abies grandis
abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries_normal.pdf",
	plot_id = "france_1001833", tree_id = 5) # 2 measurements, average growth of 7 mm/yr (realistic)

abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries_3measures_1wrong.pdf",
	plot_id = "france_725338", tree_id = 16) # 3 measurements, last is wrong

abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries_2measures_wrong.pdf",
	plot_id = "france_818095", tree_id = 4) # 2 measurements, unrealistic growth (23.8 mm/yr)

abiesGrandis = dbh_timeSeries(posterior_ls[["Abies grandis"]], plotMean = TRUE, filename = "timeSeries_germany.pdf",
	plot_id = "germany_26304_2", tree_id = 5)



example = data.table(year = rep(2000:2011, each = 3000), draw = rep(1:3000, each = 12), growth = 0)
for (i in 1:12)
	example[(3000*(i - 1) + 1):(3000*i), growth := as.numeric(qq[, , i])]

example[growth > 10, growth := NA]
example[, year := as.factor(year)]
setDF(example)

ridge_density_plot <-
	example %>%
	ggplot(aes(x = growth, y = year, fill = year)) +
	ggridges::geom_density_ridges(quantile_lines = TRUE, quantile_fun=function(x,...)mean(x), scale = 1) +
	theme(legend.position = "none", axis.title = element_text(size = 24), axis.text = element_text(size = 18),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA_character_),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("test.pdf", ridge_density_plot)


## Acer opalus
acerOpalus = dbh_timeSeries(posterior_ls[["Acer opalus"]], plotMean = TRUE, filename = "timeSeries.pdf",
	plot_id = "france_1008324", tree_id = 3, highlight_threshold_growth = 45)

acerOpalus = dbh_timeSeries(posterior_ls[["Acer opalus"]], plotMean = TRUE, filename = "timeSeries2.pdf",
	plot_id = "france_100733", tree_id = 10, highlight_threshold_growth = 45)

# Tilia platyphyllos
tiliaPlat = dbh_timeSeries(posterior_ls[["Tilia platyphyllos"]], plotMean = TRUE, filename = "timeSeries.pdf",
	plot_id = "france_121722", tree_id = 9) # index_gen = 2131, 2136

tiliaPlat = dbh_timeSeries(posterior_ls[["Tilia platyphyllos"]], plotMean = TRUE, filename = "timeSeries2.pdf",
	plot_id = "france_1004153", tree_id = 13) # index_gen = 127, 132

tiliaPlat = dbh_timeSeries(posterior_ls[["Tilia platyphyllos"]], plotMean = TRUE, filename = "timeSeries_3measures_1wrong.pdf",
	plot_id = "france_873749", tree_id = 5) # index_gen = 20330, 20335, 20340
aa = probaExtremeObs(posterior_ls[["Tilia platyphyllos"]], plot_id = "france_873749", tree_id = 5)



posteriorSim = posterior_ls[["Abies grandis"]]

aa = probaExtremeObs(posteriorSim, plot_id = "france_725338", tree_id = 16)
aa = probaExtremeObs(posteriorSim, plot_id = "france_818095", tree_id = 4)
aa = probaExtremeObs(posteriorSim, plot_id = "france_1001833", tree_id = 5)


names(posteriorSim)

draws = posteriorSim$posterior_draws
draws$metadata()$stan_variables

paste0(posteriorSim$info[["path"]], posteriorSim$info[["run"]], "_",
	c("ba_normalisation", "climate_normalisation", "dbh_normalisation", "ph_normalisation"), ".rds")

growth_fct = function(dbh, pr, tas, ph, basalArea, params, sd_dbh, standardised_params = FALSE, standardised_variables = FALSE, ...)

n_rep = as.integer(generate_quantities$info["iter_sampling"]) * as.integer(generate_quantities$info["n_chains"])

dt_dharma = data.table(rep_id = rep(1:n_obs, each = n_rep),
	rep_dbh = rep(stanData$Yobs, each = n_rep), # Observations
	latent_dbh = numeric(n_rep * n_obs), # Estimated latent dbh
	simulated_observations = numeric(n_rep * n_obs))


yearly_latent_dbh = draws$draws("yearly_latent_dbh") # iter_sampling * n_chains * n_latentStates
dim(yearly_latent_dbh)
yearly_latent_dbh_1729_1730 = sd_dbh*yearly_latent_dbh[, , indices[1729, index_gen]:indices[1730, index_gen]]

growth_array_1 = yearly_latent_dbh_1729_1730[, , 2] - yearly_latent_dbh_1729_1730[, , 1]
growth_array_2 = yearly_latent_dbh_1729_1730[, , 3] - yearly_latent_dbh_1729_1730[, , 2]
growth_array_3 = yearly_latent_dbh_1729_1730[, , 4] - yearly_latent_dbh_1729_1730[, , 3]
growth_array_4 = yearly_latent_dbh_1729_1730[, , 5] - yearly_latent_dbh_1729_1730[, , 4]
growth_array_5 = yearly_latent_dbh_1729_1730[, , 6] - yearly_latent_dbh_1729_1730[, , 5]

my_avg = mean(growth_array_1)
my_sd = sd(growth_array_1)

mu = log(my_avg^2/sqrt(my_sd^2 + my_avg^2))
sigma = sqrt(log(my_sd^2/my_avg^2 + 1))

pdf("test.pdf")
hist(growth_array_1)
curve(dlnorm(x, meanlog = mu, sdlog = sigma), lwd = 2, col = "red", add = TRUE)
dev.off()


growth_mean_logNormal = function(sigmaProc, expected_growth)
	return(log(expected_growth^2/sqrt(sigmaProc^2 + expected_growth^2)))

growth_sd_logNormal = function(sigmaProc, expected_growth)
	return (sqrt(log(sigmaProc^2/expected_growth^2 + 1)))

curve(growth_mean_logNormal(x, 9), 0, 25)


params_meanLog = growth_mean_logNormal(seq(9.997225, 12.701388, length.out = 20), 9)
params_sdLog = growth_sd_logNormal(seq(9.997225, 12.701388, length.out = 20), 9)

colours = MetBrewer::met.brewer("Hokusai3", length(params_meanLog))
colours_str = grDevices::colorRampPalette(colours)(length(params_meanLog))

curve(dlnorm(x, params_meanLog[1], params_sdLog[1]), 0, 25, col = colours_str[1], lwd = 2, xlab = "growth", ylab = "distrib",
	ylim = c(0, 0.125))

for (i in 2:length(params_meanLog))
	curve(dlnorm(x, params_meanLog[i], params_sdLog[i]), 0, 25, add = TRUE, col = colours_str[i], lwd = 2)



#* --------------------------------------------------------------------------------------------
#! ********************************************************************************************
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================    WHAT FOLLOWS ARE OLDER STUFF, SOME ARE GOOD    ================ ODOT#
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! ********************************************************************************************
#* --------------------------------------------------------------------------------------------

#! TO REMOVE AFTER !!!!!!!!!!!!!!!
results = readRDS("Abies grandis/growth-run=1-2022-05-18_00h19.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
# results = readRDS("Abies grandis/growth-run=1-2022-06-09_13h18.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
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

# draws = temporary_array
draws = yearly_latent_dbh_array
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


## Proba observation from extreme distribution
# proba_extreme_obs_array = generate_quantities$draws("proba_extreme_obs")
posteriorSim = posterior_ls[["Abies grandis"]]
proba_extreme_obs_array = posteriorSim[["posterior_draws"]]$draws("proba_extreme_obs")
range(proba_extreme_obs_array)

# With first tree
temporary_array = proba_extreme_obs_array[, , 1:2]
mean(temporary_array)

# With tree parent number 812, corresponds to indices[1729]
temporary_array = proba_extreme_obs_array[, , 1729:1730]
temporary_array1 = proba_extreme_obs_array[, , 1729]
mean(temporary_array)
temporary_array2 = proba_extreme_obs_array[, , 1730]
mean(temporary_array)

cor(temporary_array1, temporary_array2)

dd = dim(temporary_array)

colours = MetBrewer::met.brewer("Hokusai3", dd[2]*dd[3])
colours_str = grDevices::colorRampPalette(colours)(dd[2]*dd[3])

pdf("chains_proba_extreme.pdf", height = 11.25, width = 20)

plot(0, type = "n", xlim = c(0, dd[1]), ylim = c(0, 1), ylab = "proba extreme error", xlab = "Iterations")

for (param in 1:dd[3])
{
	for (chain in 1:dd[2])
		lines(1:dd[1], temporary_array[, chain, param], col = colours_str[(param - 1)*dd[2] + chain])
}
legend(x = "topleft", legend = c(paste("1729 -", 1:dd[2]), paste("1730 -", 1:dd[2])), fill = colours_str, box.lwd = 0)

dev.off()

lazyPosterior(temporary_array[,,1], filename = "1729", mean = 50, sd = 0.1, n = 10000) # Mean and sd are dummy data outside of plot window
lazyPosterior(temporary_array[,,2], filename = "1730", mean = 50, sd = 0.1, n = 10000)

mainFolder = "/home/amael/project_ssm/inventories/growth/"
treeData = readRDS(paste0(mainFolder, "standardised_european_growth_data_reshaped.rds"))
ls_species = sort(treeData[, unique(speciesName_sci)])

for (i in c(2, 5, 17, 48)) # (species in ls_species)
{
	species = ls_species[i]
	temp = treeData[speciesName_sci == species]
	temp[, next_dbh := shift(dbh, n = 1, type = "lead", fill = NA), by = .(plot_id, tree_id)]
	computeGrowth(temp)
	temp = na.omit(temp)
	temp[, dbh_75 := 0.75*dbh]
	temp = temp[next_dbh > dbh_75]
	temp = temp[growth < 20]
	# temp = temp[growth > -5]
	pdf(paste0("histG_", species, ".pdf"))
	hist(temp[, growth], breaks = seq(min(temp[, growth]), max(temp[, growth]), length.out = 20),
		main = round(x = quantile(temp[, growth], seq(0, 1, .1)), digits = 1))
	dev.off()
}

#### For a given tree, obviously measured without extreme error, is the latent average growth realistic and does it match the observation?
latG = results$draws("latent_growth")
dim(latG)
(treeData[2, dbh] - treeData[1, dbh]) # diametral increment
sum(sd_dbh*apply(latG[, , 1:5], 3, mean))
vv = numeric(3000)
for (i in 1:1000)
	for (j in 1:3)
		vv[(i - 1)*3 + j] = sum(sd_dbh*latG[i, j, 1:5])

mean(vv)
pdf("test.pdf")
hist(vv, breaks = seq(min(vv) - 1, max(vv) + 1, by = 0.5))
dev.off()



# dbh1 = c(12, 18, 7, 24, 31)
# dbh2 = c(15, 19, 17, 29, 35)
# dbh = c(dbh1, dbh2)

# years = c(5, 5, 4, 7, 2)
# gg = numeric(5)
# for (i in 1:5)
# 	gg[i] = (dbh2[i] - dbh1[i])/years[i]

# my_sd = sd(dbh)

# gg2 = numeric(5)
# for (i in 1:5)
# 	gg2[i] = (dbh2[i]/my_sd - dbh1[i]/my_sd)/years[i]

# gg/my_sd

####! CRASH TEST ZONE

## Check the proportion of unrealistic growth
# Load results
results = readRDS("Tilia platyphyllos/ruger_growth-run=1-2022-07-30_03h16_ruger.rds")
stanData = readRDS("Tilia platyphyllos/1_stanData.rds")
sd_dbh = readRDS("Tilia platyphyllos/1_dbh_normalisation.rds")[, sd]

posterior_latent_growth = sd_dbh*results$draws("latent_growth") # (iter_sampling, n_chains, n_latentGrowth)

quantile(posterior_latent_growth, seq(0, 1, 0.1))

ll = as.data.table(posterior_latent_growth) # (iter_sampling * n_chains * n_latentGrowth, 4)

threshold = 10
ll[value > threshold, .N]*100/ll[, .N]

ll[value > threshold, .N*100/3000, by = variable]

quantile(posterior_latent_growth, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1))

pdf("distrib_growth.pdf", height = 8, width = 8)
hist(ll[value > threshold, .N*100/3000, by = variable][, V1])
dev.off()

# Generate quantities for a fixed environment, dbh is the variable
gq_model = cmdstan_model("./generate_probaGrowth.stan")

# Access data
n_chains = results$num_chains()

# Change data (fixed environment, size grandient)
set.seed(1969-08-18)
indices = readRDS("Tilia platyphyllos/1_indices.rds")
patch_id = sample(x = 1:indices[.N, plot_index], size = 1)

precip = mean(stanData$precip[indices[plot_index == patch_id, unique(index_clim_start)]:
	indices[plot_index == patch_id, unique(index_clim_end)]])
temperature = mean(stanData$tas[indices[plot_index == patch_id, unique(index_clim_start)]:
	indices[plot_index == patch_id, unique(index_clim_end)]])
ph = stanData$ph[patch_id]
basalArea = mean(stanData$standBasalArea[indices[plot_index == patch_id, unique(index_clim_start)]:
	indices[plot_index == patch_id, unique(index_clim_end)]])


stanData$environment = c(precip = (precip - stanData$pr_mu)/stanData$pr_sd,
	temperature = (temperature - stanData$tas_mu)/stanData$tas_sd,
	ph = (ph - stanData$ph_mu)/stanData$ph_sd,
	basalArea = (basalArea - stanData$ba_mu)/stanData$ba_sd)

stanData$n_dbh = 1000
stanData$n_threshold = 3

stanData$dbh0 = seq(50, 550, length.out = stanData$n_dbh)
stanData$threshold = seq(10, 20, length.out = stanData$n_threshold)

# Simulations
generate_quantities = gq_model$generate_quantities(results$draws(), data = stanData, parallel_chains = n_chains)
probaGrowth = generate_quantities$draws("probaGrowth_beyondThreshold") # iter_sampling * n_chains * n_dbh * n_threshold

probaGrowth_mean = apply(probaGrowth, 3, mean)

probaGrowth_dt = data.table(dbh = rep(stanData$dbh0, stanData$n_threshold),
	threshold = rep(stanData$threshold, each = stanData$n_dbh),
	proba = probaGrowth_mean)

colours = MetBrewer::met.brewer("Hokusai3", stanData$n_threshold)
colours_str = grDevices::colorRampPalette(colours)(stanData$n_threshold)

pdf("test.pdf", height = 6, width = 17)
plot(stanData$dbh0, probaGrowth_dt[threshold == 10, proba], type = "l", col = colours_str[1], lwd = 2,
	ylim = c(0, 1.01*max(probaGrowth_mean)), xlab = "dbh", ylab = "Proba", las = 1)
lines(stanData$dbh0, probaGrowth_dt[threshold == 15, proba], col = colours_str[2], lwd = 2)
lines(stanData$dbh0, probaGrowth_dt[threshold == 20, proba], col = colours_str[3], lwd = 2)
legend(x = "topleft", legend = paste(stanData$threshold, "mm"), fill = colours_str, box.lwd = 0)
dev.off()

####! END CRASH TEST ZONE