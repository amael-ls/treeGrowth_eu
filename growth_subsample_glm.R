
#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)

#### Load data
## Common data
species = "Fagus sylvatica" # Found that run 2 and 3 disagree the most for parameter dbh^2
path = paste0("./", species, "/")
run = 4

## Data containing growth, environment, and indices
stanData = readRDS(paste0(path, run, "_stanData_classic.rds"))
treeData = readRDS(paste0(path, run, "_treeData_classic.rds"))
indices = readRDS(paste0(path, run, "_indices.rds"))[["indices_avg"]]

## For soil I need to load the whole data
soil_folder = "/home/amael/project_ssm/inventories/growth/"
if (!dir.exists(soil_folder))
	stop(paste0("Folder\n\t", soil_folder, "\ndoes not exist"))

soil = readRDS(paste0(soil_folder, "europe_reshaped_soil.rds"))
soil = soil[plot_id %in% treeData[, plot_id]]

## Remove NAs
growth_dt = na.omit(treeData)

#### Format data for GLM
## Add predictors to indices dt, rows are in the good order (see stanData in growth_subsample_classic.R)
indices[, precip := stanData$precip]
indices[, tas := stanData$tas]
indices[, standBasalArea := stanData$standBasalArea]

indices = data.table::merge.data.table(x = indices, y = soil, by = "plot_id")

if (!all.equal(indices$year, indices$year_start))
	stop("I assumed all.equal(indices$year, indices$year_start). It is not the case => need to think for merging environment and growth!")

## Remove useless columns. StandBasalArea from growth_dt is only for the measured year, while it is the avg in indices (which I want)
indices[, c("index_clim_start", "index_clim_end", "year_start", "year_end") := NULL]
growth_dt[, c("speciesName_sci", "nfi_id", "x", "y", "country", "taxonID", "deltaYear", "standBasalArea") := NULL]

## Merging data by plot and year
growth_dt = data.table::merge.data.table(x = growth_dt, y = indices, by = c("plot_id", "year"))

## Normalised and centred (except dbh that is only normalised)
growth_dt[, dbh := dbh/stanData$sd_dbh]
growth_dt[, precip := (precip - stanData$pr_mu)/stanData$pr_sd]
growth_dt[, tas := (tas - stanData$tas_mu)/stanData$tas_sd]
growth_dt[, ph := (ph - stanData$ph_mu)/stanData$ph_sd]
growth_dt[, standBasalArea := (standBasalArea - stanData$ba_mu)/stanData$ba_sd]

growth_dt[, growth := growth/stanData$sd_dbh]

## Keep only positive growths (lognormal distribution and no measurement errors)
growth_dt = growth_dt[growth > 0]

#### Run GLM
test2 = glm(growth_dt[, log(growth)] ~ 1 + growth_dt[, dbh] + growth_dt[, dbh^2] + growth_dt[, precip] + growth_dt[, precip^2] +
	growth_dt[, tas] + growth_dt[, tas^2] + growth_dt[, ph] + growth_dt[, ph^2] + growth_dt[, standBasalArea])

aa1 = summary(test1)
aa2 = summary(test2)

curve(dnorm(x, mean = aa2$coefficients["growth_dt[, dbh^2]", "Estimate"], sd = aa2$coefficients["growth_dt[, dbh^2]", "Std. Error"]),
	lwd = 2, col = "#CD1A21", from = -0.13, to = -0.08, xlab = "x", ylab = "Posterior")

curve(dnorm(x, mean = aa2$coefficients["growth_dt[, dbh^2]", "Estimate"], sd = aa2$coefficients["growth_dt[, dbh^2]", "Std. Error"]),
	lwd = 2, col = "#548912", add = TRUE)
