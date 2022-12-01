
#### Aim of prog: Create dummy data to check parameters recovery
## Comments
# The creation of these data is based on Quercus robur data. I noticed that the diameters distribution of the data is bimodal, with a
# 	clear separation between the two modes around 300 mm. Therefore, I decided to sample two gamma distributions with means and variances
#	taken from the Quercus data with dbh above or below 300 mm.

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)

#### Tool function
## Sampling N integers from 1 to X, with the constraint that the sum should be Z
# Inspired from https://stackoverflow.com/questions/3589214/generate-random-numbers-summing-to-a-predefined-value
constrained_sum_sample_pos = function(n, total)
{
	if (n > total)
		stop("The number of samples cannot be larger than the total")

	dividers = sort(sample(x = 1:(total - 1), size = n - 1, replace = FALSE))
	divTot = c(dividers, total)
	divZero = c(0, dividers)

	return (divTot - divZero)
}

#### Define variables
## Common variables
set.seed(1969-08-18) # Woodstock seed

n_indiv = 300
n_plot = 5
n_measurements = 2 # Number of measurements per individual
delta_t = 3 # Number of growing years between two measurements, e.g., 2000 --> 2003, 3 growing years 2000-2001, 2001-2002, 2002-2003

print(paste("In average there are", round(n_indiv/n_plot, 2), "individuals per plot"))

n_annual_growth_per_indiv = (n_measurements - 1)*delta_t

## Parameters (for the scaled data, so no rescaling to do)
# Intercept
beta0 = -4

# Dbh slopes
beta1 = 0.42
beta2 = -0.06 # Quadratic term

# Temperature slope
beta3 = 0.06
beta4 = -0.01 # Quadratic term

# Variance
sigmaProc = 0.42

#### Create data
## Initial diameters
treeData = data.table(dbh1 = numeric(n_indiv), big_dbh = sample(x = c(FALSE, TRUE), size = n_indiv, replace = TRUE, prob = c(0.4, 0.6)))

treeData[(big_dbh), dbh1 := rgamma(n = .N, shape = 545^2/26435, rate = 545/26435)]
treeData[!(big_dbh), dbh1 := rgamma(n = .N, shape = 196^2/3678, rate = 196/3678)]

# hist(treeData$dbh1, breaks = seq(min(treeData$dbh1) - 3, max(treeData$dbh1) + 3, 3))

treeData[, big_dbh := NULL]
sd_dbh_orig = treeData[, sd(dbh1)]

## Temperature
temperature_avg = rnorm(n = n_annual_growth_per_indiv, mean = 9.7, sd = 0.3)
temperature_sd = rgamma(n = n_annual_growth_per_indiv, shape = 0.82^2/0.0015, rate = 0.82/0.0015)

temperature = matrix(data = NA, nrow = n_plot, ncol = n_annual_growth_per_indiv)

for (year in 1:n_annual_growth_per_indiv)
	temperature[, year] = rnorm(n = n_plot, mean = temperature_avg[year], sd = temperature_sd[year])

mu_temp = mean(temperature)
sd_temp = sd(temperature)

temperature = as.data.table(temperature)
setnames(temperature, paste0("temperature", 1:n_annual_growth_per_indiv))
temperature[, plot_id := 1:n_plot]

## Assign trees to plots
nb_trees_per_plot = constrained_sum_sample_pos(n = n_plot, total = n_indiv)
treeData[, plot_id := rep(1:n_plot, times = nb_trees_per_plot)]
treeData[, tree_id := 1:.N, by = plot_id]

all.equal(treeData[, .N, by = plot_id][, N], nb_trees_per_plot)

## Merge treeData and temperatures
treeData = treeData[temperature, on = "plot_id"]

## Generate all the diameters
for (year in 1:n_annual_growth_per_indiv)
{
	current_temperature = paste0("temperature", year)
	current_dbh = paste0("dbh", year)
	next_dbh = paste0("dbh", year + 1)

	log_reg = beta0 + beta1*treeData[[current_dbh]]/sd_dbh_orig + beta2*(treeData[[current_dbh]]/sd_dbh_orig)^2 +
		beta3*(treeData[[current_temperature]] - mu_temp)/sd_temp + beta4*(treeData[[current_temperature]] - mu_temp)^2/sd_temp^2
	current_growth = rlnorm(n = n_indiv, meanlog = log_reg + log(sd_dbh_orig), sdlog = sigmaProc)

	treeData[, c(next_dbh) := treeData[[current_dbh]] + current_growth]

	print(range(current_growth))
}

## Generate the observed growth (recorded only at frequency delta_t)
dbh_recorded = paste0("dbh", seq(1, n_annual_growth_per_indiv + 1, by = delta_t))

for (i in seq_along(dbh_recorded[-length(dbh_recorded)]))
{
	current_dbh = dbh_recorded[i]
	next_rec_dbh = dbh_recorded[i + 1]
	obs_growth = paste0("growth", i)

	treeData[, c(obs_growth) := (treeData[[next_rec_dbh]] - treeData[[current_dbh]])/delta_t]
}


setcolorder(treeData, c("plot_id", "tree_id", paste0("dbh", 1:(n_annual_growth_per_indiv + 1)),
	paste0("growth", seq_along(dbh_recorded[-length(dbh_recorded)])),
	paste0("temperature", 1:n_annual_growth_per_indiv)))

saveRDS(list(
		treeData = treeData,
		parameters = c(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4, sigmaProc = sigmaProc),
		scaling = c(sd_dbh_orig = sd_dbh_orig, mu_temp = mu_temp, sd_temp = sd_temp),
		infos = c(n_indiv = n_indiv, n_plot = n_plot, n_measurements = n_measurements, delta_t = delta_t,
			n_annual_growth_per_indiv = n_annual_growth_per_indiv)
	), "dummyData.rds")
