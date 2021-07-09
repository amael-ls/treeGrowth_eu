
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



###########################################################?
######## 		Second PART: Create indices 		########
###########################################################?

#### Tool function
fillYears = function(years)
{
	if (length(years) < 2)
		stop("From fillYears: Their should be at least two years to fill the gaps")

	if (is.unsorted(years))
		stop("From fillYears: years are assumed to be sorted!")

	fill_years = years[1]:years[length(years)]
	indices = which(fill_years %in% years)
	
	return (list(fill_years = fill_years, indices = indices))
}

#### Get precipitations (in real, from a raster but here from treeStates_dt)
values = unique(treeStates_dt[, .(plot_id, year, precipitations)])
setnames(values, old = "plot_id", new = "id")

#### 'Joining' climate with tree data
count = 0
start = 0
end = 0
iter = 0

nbIndiv = unique(treeData[, .(plot_id, tree_id)])[, .N]
length_filled_years = sum(treeData[, max(year) - min(year) + 1, by = .(plot_id, tree_id)][, V1])

indices = data.table(year = integer(length_filled_years), tree_id = integer(length_filled_years),
	plot_id = integer(length_filled_years), index_gen = integer(length_filled_years),
	index_precip_start = integer(length_filled_years), index_precip_end = integer(length_filled_years))

for (plot in treeData[, unique(plot_id)])
{
	for (indiv in treeData[plot_id == plot, unique(tree_id)])
	{
		years_indices = fillYears(treeData[plot_id == plot & tree_id == indiv, year])
		start = end + 1
		end = end + length(years_indices[["fill_years"]])
		indices[start:end, year := years_indices[["fill_years"]]]
		indices[start:end, tree_id := indiv]
		indices[start:end, plot_id := plot]
		indices[years_indices[["indices"]] + count, index_gen := years_indices[["indices"]] + count]
		count = count + years_indices[["indices"]][length(years_indices[["indices"]])]
		iter = iter + 1
		if (iter %% 1000 == 0)
			print(paste(round(iter*100/nbIndiv, digits = 3), "% done"))
	}
}

## Create the indices for the climate, and format climate data
start = 0
end = 0
count = 0
iter = 0

length_clim = sum(treeData[, max(year) - min(year) + 1, by = .(plot_id)][, V1])
precipitations_yearly = numeric(length_clim)

for (plot in indices[, unique(plot_id)])
{
	precip_years = indices[plot_id == plot, sort(unique(year))]
	
	start = end + 1
	end = start + length(precip_years) - 1

	precipitations_yearly[start:end] = as.numeric(values[id == plot & year %in% precip_years , precipitations])

	for (tree in indices[plot_id == plot, unique(tree_id)])
	{
		precip_start = min(which(indices[tree_id == tree & plot_id == plot, year] %in% precip_years)) + count
		precip_end = max(which(indices[tree_id == tree & plot_id == plot, year] %in% precip_years)) + count

		indices[tree_id == tree & plot_id == plot,
			c("index_precip_start", "index_precip_end") := .(precip_start, precip_end)]
		
		iter = iter + 1
		if (iter %% 1000 == 0)
			print(paste(round(iter*100/nbIndiv, digits = 3), "% done"))
	}
	count = count + length(precip_years)
}

indices = indices[index_gen != 0]
