
# This prog generates data for toy.stan, toyPara_GPUs.stan and toyPara_reduce_sum.stan
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var

####
rm(list = ls())
graphics.off()

library(data.table)

# library(stringi)
# library(raster)
# library(sp)

options(max.print = 500)
###########################################################?
######## 		FIRST PART: Create the data 		########
###########################################################?

#### Variables
## Comon variables
# Set seed
set.seed(1969-08-18) # Woodstock seed

# Number of plots
n_plots = 4

# Number of census per plot
record_min = 4
record_max = 8

n_census_per_plot = sample(x = record_min:record_max, size = n_plots, replace = TRUE)

# Number of individuals per plot
indiv_min = 10
indiv_max = 30
n_indiv_per_plot = sample(x = indiv_min:indiv_max, size = n_plots, replace = TRUE)

# Number of data
nb_measures = sum(n_indiv_per_plot * n_census_per_plot)

## Storing data
treeStates_dt = data.table(year = integer(nb_measures), tree_id = integer(nb_measures),
	plot_id = integer(nb_measures), true_dbh = numeric(nb_measures), precipitations = numeric(nb_measures))

## Parameters
intercept = 2.4

slope_precip = 0.004
quad_slope_precip = -0.00001

slope_dbh = 1.1

sigma_process = 4
sigma_measure = 1

#### Generate complete data
init_year_plot = sample(1990:2005, size = n_plots, replace = TRUE)

count = 0
for (xy in 1:n_plots)
{
	census_years = rep(init_year_plot[xy]:(init_year_plot[xy] + n_census_per_plot[xy] - 1), n_indiv_per_plot[xy])
	nb_measures_plot = n_census_per_plot[xy] * n_indiv_per_plot[xy]
	treeStates_dt[(count + 1):(count + nb_measures_plot), c("year", "plot_id", "tree_id") :=
		.(census_years, xy, rep(1:n_indiv_per_plot[xy], each = n_census_per_plot[xy]))]

	precip = runif(n = n_census_per_plot[xy], min = 650, max = 1200)

	treeStates_dt[(count + 1):(count + nb_measures_plot), precipitations := rep(precip, n_indiv_per_plot[xy])]

	# Initial state
	initial_dbh = rgamma(n_indiv_per_plot[xy], shape = 150^2/300, rate = 150/300)
	index = seq(from = 1, to = nb_measures_plot, by = n_census_per_plot[xy]) + count
	treeStates_dt[index, true_dbh := initial_dbh]

	for (j in 1:(n_census_per_plot[xy] - 1))
	{
		next_dbh = intercept + slope_dbh*treeStates_dt[index + j - 1, true_dbh] +
			slope_precip*treeStates_dt[index + j - 1, precipitations] +
			quad_slope_precip*(treeStates_dt[index + j - 1, precipitations])^2
		treeStates_dt[index + j, true_dbh := rnorm(n = n_indiv_per_plot[xy], mean = next_dbh, sd = sigma_process)]
	}
	
	count = count + nb_measures_plot
}

treeStates_dt[, range(true_dbh)]

treeStates_dt[, observed_dbh := rnorm(n = .N, mean = true_dbh, sd = sigma_measure)]

#### Generate partial data set (keep only first and last measurements)
kept_rows = sort(c(treeStates_dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1],
	treeStates_dt[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]))
treeData = treeStates_dt[kept_rows]
setorder(treeData, plot_id, tree_id, year)



