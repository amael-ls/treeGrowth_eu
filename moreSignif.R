
#### Aim of prog: To determine, for the quadratic term, if they are 'more significant' with SSM
## Comment
#	I expect the quadratic terms to be more negative with SSM compare to classic, i.e., 'more significant'
#	Note that I only talk about negative because this is also the sign I expect for the quadratic terms


library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
source("./toolFunctions.R")

#### List species
ls_species = c("Abies alba", "Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus sylvestris", "Quercus petraea")
currentSpecies = ls_species[1]

quadratic_terms = c("dbh_slope2", "pr_slope2", "tas_slope2", "ph_slope2")

run = 1

dt = data.table(species = rep(ls_species, each = length(quadratic_terms)), param = rep(quadratic_terms, length(ls_species)))
setkey(dt, species, param)

ssm_cols = c("ssm_05", "ssm_50", "ssm_95")
classic_cols = c("classic_05", "classic_50", "classic_95")

for (currentSpecies in ls_species[1:5])
{
	tree_path = paste0("./", currentSpecies, "/")

	info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
	ssm = readRDS(paste0(tree_path, info_lastRun[["file"]]))

	info_lastRun = getLastRun(path = tree_path, begin = "^growth-", extension = "_classic.rds$", format = "ymd", run = run, hour = TRUE)
	classic = readRDS(paste0(tree_path, info_lastRun[["file"]]))

	for (quadTerm in quadratic_terms)
	{
		values_ssm = ssm$draws(quadTerm)
		values_classic = classic$draws(quadTerm)

		dt[.(currentSpecies, quadTerm), c(ssm_cols) := as.list(quantile(values_ssm, probs = c(0.05, 0.5, 0.95)))]
		dt[.(currentSpecies, quadTerm), ssm_mean := mean(values_ssm)]
		dt[.(currentSpecies, quadTerm), c(classic_cols) := as.list(quantile(values_classic, probs = c(0.05, 0.5, 0.95)))]
		dt[.(currentSpecies, quadTerm), classic_mean := mean(values_classic)]

	}
}

dt["Fagus sylvatica", .(ssm_mean, classic_mean)]

dt[(ssm_50 < 0) & (classic_50 < 0), sum(ssm_50 < classic_50)]/dt[(ssm_50 < 0) & (classic_50 < 0), .N]*100 # 61.53846
dt[(ssm_50 < classic_50) & (classic_50 < 0), .N]/dt[(ssm_50 < 0) & (classic_50 < 0), .N]*100 # 61.53846
dt[ssm_50 < 0, .N]/dt[classic_50 < 0, .N] # 15/14 = 1.071429
dt[ssm_95 < 0, .N]/dt[classic_95 < 0, .N] # 13/10 = 1.3
dt[ssm_05 > 0, .N]/dt[classic_05 > 0, .N] # 2/4 = 0.5
