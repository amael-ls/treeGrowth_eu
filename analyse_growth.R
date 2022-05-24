
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

#! CHANGE SIGMA PROC ?????

# params_names = paste0("latent_dbh_parents[", 812:813, "]")

# params_names = results$metadata()$stan_variables
# params_names = params_names[params_names != "averageGrowth"]

# for (i in 1:length(rm_names))
# 	params_names = params_names[!stri_detect(params_names, regex = rm_names[i])]

# draws_array = 143.1517*results$draws(params_names)
# dimnames(draws_array)$variable = 812:813
# np = nuts_params(results)
# str(np)
# levels(np$Parameter)

# color_scheme_set("darkgray")
# div_style = parcoord_style_np(div_color = "green", div_size = 0.05, div_alpha = 0.4)
# pdf("test.pdf", height = 11.25, width = 20)
# bayesplot::mcmc_parcoord(draws_array, size = 0.25, alpha = 0.1, np = np, np_style = div_style)
# dev.off()

# low_rhat = c("latent_dbh_parents[812]", "latent_dbh_parents[829]", "latent_dbh_parents[861]", "latent_dbh_parents[862]")

#  r$> parents_index[812]
#  [1] 1729
#  
#  r$> indices[1729:1730]
#      year tree_id       plot_id index_gen index_clim_start index_clim_end plot_index    nfi nfi_index   type nbYearsGrowth
#     <int>   <int>        <char>     <int>            <int>          <int>      <int> <char>     <int> <char>         <int>
#? 1:  2008       4 france_818095      5397           525235         525245        188 france         1 parent             5
#  2:  2013       4 france_818095      5402           525235         525245        188 france         1  child             5
#  
#  r$> treeData[1729:1730]
#     speciesName_sci nfi_id       plot_id tree_id  year      dbh        x        y standBasalArea country        taxonID
#              <char> <char>        <char>   <int> <int>    <num>    <num>    <num>          <num>  <char>         <char>
#? 1:   Abies grandis FR IFN france_818095       4  2008 511.2536 1.876306 45.76976       46.76220  france wfo-0000511178
#  2:   Abies grandis FR IFN france_818095       4  2013 630.2536 1.876306 45.76976       54.71816  france wfo-0000511178

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

params_dt = data.table(parameters = c("averageGrowth_mu", "averageGrowth_sd", "dbh_slope", "pr_slope", "pr_slope2",
	"tas_slope", "tas_slope2", "ph_slope", "ph_slope2", "competition_slope"),
	priors = c(dnorm, dgamma, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm, dnorm),
	arg1 = c(-4, 1/100, 0, 0, 0, 0, 0, 0, 0, 0),
	arg2 = c(10, 1/100, 5, 5, 5, 5, 5, 5, 5, 5),
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
}

error_dt = rbindlist(error_ls, idcol = "speciesName_sci")
correl_dt = rbindlist(correl_ls, idcol = "speciesName_sci")

saveRDS(error_dt, "./error_species.rds")
saveRDS(correl_dt, "./correlation_energy_species.rds")

plot_correl_error(error_dt, correl_dt, threshold_correl = 0.2, rm_correl = "lp__")


#! TO REMOVE AFTER !!!!!!!!!!!!!!!
results = readRDS("Abies grandis/growth-run=1-2022-05-18_00h19.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
stanData = readRDS("Abies grandis/1_stanData.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
indices = readRDS("Abies grandis/1_indices.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
treeData = readRDS("Abies grandis/1_treeData.rds") #! TO REMOVE AFTER !!!!!!!!!!!!!!!
divergences = which(results$sampler_diagnostics()[, , "divergent__"] == 1) #! TO REMOVE AFTER !!!!!!!!!!!!!!!
sd_dbh = readRDS("Abies grandis/1_dbh_normalisation.rds")[, sd] #! TO REMOVE AFTER !!!!!!!!!!!!!!!

#* --------------------------------------------------------------------------------------------
#! ********************************************************************************************
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TODO ================    WHAT FOLLOWS ARE OLDER STUFF, SOME ARE GOOD    ================ ODOT#
#? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! ********************************************************************************************
#* --------------------------------------------------------------------------------------------


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
data = treeData[1729:1730]

ll = dbhTrajectories(draws, data, plotMean = TRUE, divergences = divergences, filename = "test.pdf", highlight_max_growth = TRUE,
	highlight_min_start = TRUE, highlight_max_start = TRUE, highlight_min_end = TRUE, highlight_max_end = TRUE,
	highlight_threshold_growth = 25)

iter_mag_growth = ll[["highlight_max_growth"]]["iter_id"]

params = c("dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2", "competition_slope",
	"averageGrowth_mu", "averageGrowth_sd", "sigmaProc", "sigmaObs", "etaObs", "proba")

aa = results$draws(params)[664:668,,]

for (param in params)
	lazyTrace(draws = results$draws(param), filename = paste0("caca_", param, ".pdf"), iter1 = iter_mag_growth)

gq_draws = generate_quantities$draws("latentG_residuals")
re_draws = results$draws("latent_growth")

expectedG = sd_dbh*(gq_draws + re_draws)
qq = expectedG[664:668,,indices[1729, index_gen]:indices[1730, index_gen]]

qq[,,1:2]
qq[,,4:6]

max_state = which.max(expectedG)
local_iter = max_state %% (dim(expectedG)[1] * dim(expectedG)[2])
state = (max_state - local_iter)/(dim(expectedG)[1] * dim(expectedG)[2]) + 1

expectedG[,,state]

# r$> indices[1599:1601]
#     year tree_id       plot_id index_gen index_clim_start index_clim_end plot_index    nfi nfi_index   type nbYearsGrowth
#    <int>   <int>        <char>     <int>            <int>          <int>      <int> <char>     <int> <char>         <int>
# 1:  2007       1 france_738952      4979           474738         474748        181 france         1 parent            10
# 2:  2012       1 france_738952      4984           474738         474748        181 france         1  child            10
# 3:  2017       1 france_738952      4989           474738         474748        181 france         1  child            10

# r$> treeData[(plot_id == "france_738952") & (tree_id == 1)]
#    speciesName_sci nfi_id       plot_id tree_id  year      dbh        x        y standBasalArea country        taxonID
#             <char> <char>        <char>   <int> <int>    <num>    <num>    <num>          <num>  <char>         <char>
# 1:   Abies grandis FR IFN france_738952       1  2007 284.2423 4.429703 44.95929       35.44518  france wfo-0000511178
# 2:   Abies grandis FR IFN france_738952       1  2012 331.0423 4.429703 44.95929       48.10560  france wfo-0000511178
# 3:   Abies grandis FR IFN france_738952       1  2017 401.7071 4.429703 44.95929       30.44573  france wfo-0000511178