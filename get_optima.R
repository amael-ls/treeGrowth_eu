
#### Aim of prog: Rewrite growth as a Gaussian (when possible) and get the optimal growth conditions and sd around the optimum
## Comment:
# When a quadratic term is negative, the growth function can be rewritten as a Gaussian. This has the advantage of reparametrising the
#	growth as a function of 'mean' (i.e., the optimal explanatory value) and 'sd' around the explanatory value. Note that the terms mean and sd
#	are not the most adapted terms here, as they should be used for random variables only. Hence, I use the following wording:
#		1. The mean is the optimum
#		2. The sd is the tolerance
#
# There are more details in the appendix of the associated article
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
## Source functions
source("toolFunctions.R")

#### Load results
## Common variables
args = commandArgs(trailingOnly = TRUE)
for (i in seq_along(args))
	print(paste0("Arg ", i, ": <", args[i], ">"))

if (length(args) != 2)
	stop("Supply in this order the species and the run_id", call. = FALSE)

(species = as.character(args[1]))
(run = as.integer(args[2]))

tree_path = paste0("./", species, "/")
if (!dir.exists(tree_path))
	stop(paste0("Path not found for species <", species, ">."))

params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2", "ph_slope", "ph_slope2",
	"competition_slope") # Note that competition_slope is useless, but I still ask for it

## Results
info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
classic = readRDS(paste0(tree_path, info_lastRun[["file"]]))

ssm_params = getParams(model_cmdstan = ssm, params_names = params, type = "all")
classic_params = getParams(model_cmdstan = classic, params_names = params, type = "all")

## Load scalings
# Load climate and ph scalings (dbh and basal area useless here, as BA already standardised, and the good sd_dbh is in stanData)
sd_dbh_ssm = readRDS(paste0(tree_path, run, "_dbh_normalisation.rds"))
climate_scaling_ssm = readRDS(paste0(tree_path, run, "_climate_normalisation.rds"))
ph_scaling_ssm = readRDS(paste0(tree_path, run, "_ph_normalisation.rds"))

sd_dbh_classic = readRDS(paste0(tree_path, run, "_dbh_normalisation_classic.rds"))
climate_scaling_classic = readRDS(paste0(tree_path, run, "_climate_normalisation_classic.rds"))
ph_scaling_classic = readRDS(paste0(tree_path, run, "_ph_normalisation_classic.rds"))

scaling_ssm = rbindlist(list(sd_dbh_ssm, climate_scaling_ssm, ph_scaling_ssm))
scaling_classic = rbindlist(list(sd_dbh_classic, climate_scaling_classic, ph_scaling_classic))

scaling_dt = rbindlist(list(ssm = scaling_ssm, classic = scaling_classic), idcol = "type")
scaling_dt[, variable := stri_replace_all(str = variable, replacement = "", regex = "_avg$")] # Change the names for ease of usage

setkey(scaling_dt, type, variable)

#### Get the mean and sd for each combination of parameters
## Dimensions (according to SSM, but that does not matter: Classic were run with the same number of iterations and chains)
n_chains = ssm$num_chains()
iter_sampling = ssm$metadata()$iter_sampling

gauss_ssm = data.table(dbh = rep(-Inf, n_chains*iter_sampling), pr = -Inf, tas = -Inf, ph = -Inf)
gauss_classic = data.table(dbh = rep(-Inf, n_chains*iter_sampling), pr = -Inf, tas = -Inf, ph = -Inf)
counter = 0

for (cc in seq_len(n_chains))
{
	for (i in seq_len(iter_sampling))
	{
		ssm_temp = as.numeric(ssm_params[i, cc, ])
		names(ssm_temp) = params
		classic_temp = as.numeric(classic_params[i, cc, ])
		names(classic_temp) = params

		counter = counter + 1
		gauss_ssm[counter, c("dbh", "pr", "tas", "ph") := as.list(gaussStyle(params = ssm_temp)[["mu"]])]
		gauss_classic[counter, c("dbh", "pr", "tas", "ph") := as.list(gaussStyle(params = classic_temp)[["mu"]])]
	}
}

## Rescale the optima
gauss_ssm[, dbh := scaling_dt[.("ssm", "dbh"), sd]*dbh]
gauss_ssm[, pr := scaling_dt[.("ssm", "pr"), sd]*pr + scaling_dt[.("ssm", "pr"), mu]]
gauss_ssm[, tas := scaling_dt[.("ssm", "tas"), sd]*tas + scaling_dt[.("ssm", "tas"), mu]]
gauss_ssm[, ph := scaling_dt[.("ssm", "ph"), sd]*ph + scaling_dt[.("ssm", "ph"), mu]]

gauss_classic[, dbh := scaling_dt[.("classic", "dbh"), sd]*dbh]
gauss_classic[, pr := scaling_dt[.("classic", "pr"), sd]*pr + scaling_dt[.("classic", "pr"), mu]]
gauss_classic[, tas := scaling_dt[.("classic", "tas"), sd]*tas + scaling_dt[.("classic", "tas"), mu]]
gauss_classic[, ph := scaling_dt[.("classic", "ph"), sd]*ph + scaling_dt[.("classic", "ph"), mu]]

#### Plot the optima with the 5 - 95% interval
## Get the quantile intervals
qt_ssm = data.table(variable = colnames(gauss_ssm), p05 = -Inf, med = -Inf, p95 = -Inf, key = "variable")
qt_classic = data.table(variable = colnames(gauss_classic), p05 = -Inf, med = -Inf, p95 = -Inf, key = "variable")

for (current_col in colnames(gauss_ssm))
{
	n_NA_ssm = sum(is.na(gauss_ssm[, ..current_col]))
	n_NA_classic = sum(is.na(gauss_classic[, ..current_col]))
	if (n_NA_ssm > 0)
		warning(paste0(round(n_NA*100/gauss_ssm[, .N], 2), "% of the data for SSM were ignored because they were NA"))
	
	if (n_NA_classic > 0)
		warning(paste0(round(n_NA*100/gauss_ssm[, .N], 2), "% of the data for Classic were ignored because they were NA"))
	
	qt_ssm[current_col, c("p05", "med", "p95") := as.list(quantile(gauss_ssm[, ..current_col], c(0.05, 0.5, 0.95), na.rm = TRUE))]
	qt_classic[current_col, c("p05", "med", "p95") := as.list(quantile(gauss_classic[, ..current_col], c(0.05, 0.5, 0.95), na.rm = TRUE))]
}

saveRDS(qt_ssm, paste0(tree_path, "optima_ssm.rds"))
saveRDS(qt_classic, paste0(tree_path, "optima_classic.rds"))
