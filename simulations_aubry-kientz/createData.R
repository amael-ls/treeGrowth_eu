
#### Aim of prog: Create dummy data to check parameters recovery
## Comments
# The creation of these data is based on Quercus robur data. I noticed that the diameters distribution of the data is bimodal, with a
# 	clear separation between the two modes around 300 mm. Therefore, I decided to sample two gamma distributions with means and variances
#	taken from the Quercus data with dbh above or below 300 mm.
# You need to run this program with n_measurements = 2, 3, 4, 5, 7, and 13, respectively
#
# 19th June 2024: I added a parameter for correlated climate data within plot or not. So far, it was correlated. Now, I want to see what happens
#	with non correlated plots and less climate diverstity

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(stringi)

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

	return(divTot - divZero)
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
	abs(x - round(x)) < tol

#### Define variables
## Common variables
set.seed(1969 - 08 - 18) # Woodstock seed

n_indiv = 400
n_plot = 50 # round(n_indiv/5) # 20 individuals per plot in average
n_measurements = 2 # Number of dbh measurements per individual. I chose among (2:5, 7, 13), to generate 1, 2, 3, 4, 6 and 12 growth observations
length_time_series = 12 # For how many years
delta_t = length_time_series/(n_measurements - 1) # Frequence of measurements

print(paste("delta_t =", delta_t))

if (!is.wholenumber(delta_t))
	stop("delta_t is not an integer... It would make more sense it is, as NFI data are measured with an integer frequency")

print(paste("In average there are", round(n_indiv/n_plot, 2), "individuals per plot"))

n_annual_growth_per_indiv = (n_measurements - 1)*delta_t # Should be length_time_series by definition

## Parameters (for the scaled data, so no rescaling to do)
# Intercept
beta0 = -4

# Dbh slopes
beta1 = 0.42
beta2 = -0.06 # Quadratic term

# Temperature slope
beta3 = 0.076
beta4 = -0.01 # Quadratic term

# Variance
sigmaProc = 0.42
unscaled_sigmaObs = 0.1

#### Create data
## Initial diameters
treeData = data.table(dbh1 = numeric(n_indiv), big_dbh = sample(x = c(FALSE, TRUE), size = n_indiv, replace = TRUE, prob = c(0.4, 0.6)))

treeData[(big_dbh), dbh1 := rgamma(n = .N, shape = 545^2/26435, rate = 545/26435)]
treeData[!(big_dbh), dbh1 := rgamma(n = .N, shape = 196^2/3678, rate = 196/3678)]

if (!file.exists("./histogram_dbh.tex"))
{
	tikz("./histogram_dbh.tex", height = 3, width = 3)
	hist(treeData$dbh1, breaks = seq(min(treeData$dbh1) - 30, max(treeData$dbh1) + 30, 30), xlab = "Diameter", ylab = "Frequency",
		main = "", las = 1)
	dev.off()
}

treeData[, big_dbh := NULL]
sd_dbh_orig = treeData[, sd(dbh1)]

sigmaObs = unscaled_sigmaObs/sd_dbh_orig

## Temperature
temperature_avg = rnorm(n = n_plot, mean = 9, sd = 1) # Average temperature of each plot
print(paste("Range temperature avg:", paste(round(range(temperature_avg), 2), collapse = " - ")))
temperature_sd = rgamma(n = n_plot, shape = 0.5^2/0.015, rate = 0.5/0.015) # sd for year to year variation for each plot
print(paste("Range temperature sd:", paste(round(range(temperature_sd), 3), collapse = " - ")))

temperature = matrix(data = NA, nrow = n_plot, ncol = n_annual_growth_per_indiv) # n_plot series of temperature, of length n_annual_growth_indiv

correlated = FALSE

if (correlated)
{
	for (year in 1:n_annual_growth_per_indiv)
		temperature[, year] = rnorm(n = n_plot, mean = temperature_avg, sd = temperature_sd)
} else {
	temperature = matrix(data = rnorm(n = n_plot*n_annual_growth_per_indiv, mean = mean(temperature_avg), sd = max(temperature_sd)),
		nrow = n_plot, ncol = n_annual_growth_per_indiv)
}


mu_temp = mean(temperature)
sd_temp = sd(temperature)

temperature = as.data.table(temperature)
setnames(temperature, paste0("temperature", 1:n_annual_growth_per_indiv))
temperature[, plot_id := 1:n_plot]

# Get the columns that match the pattern
temperatur_cols = grep("temperature", names(temperature), value = TRUE)

# Apply the function to each row of the selected columns
temperature[, avg := apply(.SD, 1, mean), .SDcols = temperatur_cols]

## Assign trees to plots
nb_trees_per_plot = constrained_sum_sample_pos(n = n_plot, total = n_indiv)
treeData[, plot_id := rep(1:n_plot, times = nb_trees_per_plot)]
treeData[, tree_id := seq_len(.N), by = plot_id]

all.equal(treeData[, .N, by = plot_id][, N], nb_trees_per_plot)

if (correlated && !file.exists("./histogram_nbPerPlot.tex"))
{
	tikz("./histogram_nbPerPlot.tex", height = 3, width = 3)
	hist(nb_trees_per_plot, breaks = seq(0, max(nb_trees_per_plot), 1), xlab = "Number of indiviudals per plot",
		ylab = "Frequency", main = "", las = 1)
	dev.off()
}

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

for (i in seq_along(dbh_recorded))
{
	current_dbh = dbh_recorded[i]
	obs_dbh = paste0("obs_", current_dbh)

	treeData[, c(obs_dbh) := rnorm(n = .N, mean = treeData[[current_dbh]], sd = unscaled_sigmaObs)]
}

for (i in seq_along(dbh_recorded[-length(dbh_recorded)]))
{
	current_dbh = dbh_recorded[i]
	next_rec_dbh = dbh_recorded[i + 1]

	obs_current_dbh = paste0("obs_", current_dbh)
	obs_next_rec_dbh = paste0("obs_", next_rec_dbh)

	obs_growth = paste0("obs_growth", i)
	avg_temperature = paste0("avg_temperature", i)

	id_dbh_current = as.integer(stri_sub(str = current_dbh, from = stri_locate_first(str = current_dbh, regex = "[:digit:]")[, "end"]))
	id_dbh_next = as.integer(stri_sub(str = next_rec_dbh, from = stri_locate_first(str = next_rec_dbh, regex = "[:digit:]")[, "end"])) - 1
	temperature_cols = paste0("temperature", id_dbh_current:id_dbh_next)

	treeData[, c(obs_growth) := (treeData[[obs_next_rec_dbh]] - treeData[[obs_current_dbh]])/delta_t]
	treeData[, c(avg_temperature) := sum(.SD)/(delta_t*.N), .SDcols = temperature_cols, by = .(plot_id)]
	# The *.N is to compensate that I do not compute per .(plot_id, tree_id)! I divide by the number of tree_id
}

setcolorder(treeData, c("plot_id", "tree_id", paste0("dbh", 1:(n_annual_growth_per_indiv + 1)),
	paste0("obs_", dbh_recorded),
	paste0("obs_growth", seq_along(dbh_recorded[-length(dbh_recorded)])),
	paste0("temperature", 1:n_annual_growth_per_indiv)))

if (correlated)
{
	dataFilename = paste0("dummyData_plot=", n_plot, "_indiv=", n_indiv, "_deltaT=", delta_t, ".rds")
} else {
	dataFilename = paste0("dummyData_plot=", n_plot, "_indiv=", n_indiv, "_deltaT=", delta_t, "_nonCorrelated.rds")
}

saveRDS(list(
		treeData = treeData,
		parameters = c(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4, sigmaProc = sigmaProc),
		scaling = c(sd_dbh_orig = sd_dbh_orig, mu_temp = mu_temp, sd_temp = sd_temp),
		infos = c(n_indiv = n_indiv, n_plot = n_plot, n_measurements = n_measurements, delta_t = delta_t,
			n_annual_growth_per_indiv = n_annual_growth_per_indiv),
		sigmaObs = sigmaObs
	), dataFilename)

print(paste0("Data saved under the name <", dataFilename, ">"))
