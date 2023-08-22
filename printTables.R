
#### Aim of prog: Print tables, useful to create tables in latex (appendix app_S7.tex)

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)

#### Tool functions
source("toolFunctions.R")

#### Print tables
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")

params_names = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope")

params_values = printParams(ls_species = ls_species, params_names = params_names, run = 1)

ssm_05 = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("q5"))
ssm_med = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("med"))
ssm_95 = data.table::dcast(params_values[["ssm_params"]], parameter ~ speciesName_sci, value.var = c("q95"))

classic_05 = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("q5"))
classic_med = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("med"))
classic_95 = data.table::dcast(params_values[["classic_params"]], parameter ~ speciesName_sci, value.var = c("q95"))

data.table:::print.data.table(x = ssm_med, digits = 2)
data.table:::print.data.table(x = classic_med, digits = 2)
