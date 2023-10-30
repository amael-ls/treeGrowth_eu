
#### Aim of prog: Compute and plot predicted vs observed

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)

#### Tool functions
source("toolFunctions.R")

## Indexing for predict
indices_predict = function(treeData, climate, time_space)
{
	#### Tool functions
	## Function to fill the gap between two years and get the index of the provided years
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

	## Function to create the individual indices for the state-space model (stan)
	indices_state_space = function(trees_NFI)
	{
		count = 0
		start = 0
		end = 0
		iter = 0

		if (length(trees_NFI[, unique(speciesName_sci)]) > 1)
			stop(paste0("from indices_state_space: Their should be only one species of trees. Currently there are: \n- ",
				paste0(trees_NFI[, unique(speciesName_sci)], collapse = "\n- ")))
		
		nbIndiv = unique(trees_NFI[, .(plot_id, tree_id)])[, .N]
		length_filled_years = sum(trees_NFI[, max(year) - min(year) + 1, by = .(plot_id, tree_id)][, V1])

		indices = data.table(year = integer(trees_NFI[, .N]), tree_id = character(trees_NFI[, .N]),
			plot_id = character(trees_NFI[, .N]), index_gen = integer(trees_NFI[, .N]),
			index_clim_start = integer(trees_NFI[, .N]), index_clim_end = integer(trees_NFI[, .N]))

		for (plot in trees_NFI[, unique(plot_id)])
		{
			for (indiv in trees_NFI[plot, unique(tree_id)])
			{
				years_indices = fillYears(trees_NFI[.(plot, indiv), year])
				
				start = end + 1
				end = start + length(years_indices[["indices"]]) - 1
				
				indices[start:end, year := years_indices[["fill_years"]][years_indices[["indices"]]]]
				indices[start:end, tree_id := indiv]
				indices[start:end, plot_id := plot]
				indices[start:end, index_gen := years_indices[["indices"]] + count]

				count = count + years_indices[["indices"]][length(years_indices[["indices"]])]
				iter = iter + 1
				if (iter %% 1000 == 0)
					print(paste(round(iter*100/nbIndiv, digits = 1), "% done"))
			}
		}
		print("100% done")
		return (indices)
	}

	## Function to create the individual indices to match with climate (stan)
	# This function will modify the indices data table by reference and avoiding a copy.
	indices_climate = function(indices_dt, climate, time_space)
	{
		for (plot in indices_dt[, unique(plot_id)])
		{
			min_year = time_space[plot, min_year]
			max_year = time_space[plot, max_year]
			tree = indices_dt[plot, unique(tree_id)][1]
			for (tree in indices_dt[plot, unique(tree_id)])
			{
				clim_start = climate[.(plot, min_year), row_id]
				clim_end = climate[.(plot, max_year), row_id]

				indices_dt[.(plot, tree),
					c("index_clim_start", "index_clim_end") := .(clim_start, clim_end)]
			}
		}
		print("100% done")
	}

	## Function to create the (plot-year) indices to compute average climate
	# This function will modify the indices data table by reference and avoiding a copy.
	indices_climate_avg = function(indices_dt, climate)
	{
		for (plot in indices_dt[, unique(plot_id)])
		{
			all_years = indices_dt[plot, unique(year)]

			for (i in 1:(length(all_years) - 1))
			{
				clim_start = climate[.(plot, all_years[i]), row_id]
				clim_end = climate[.(plot, all_years[i + 1]), row_id]

				indices_dt[.(plot, all_years[i]),
					c("index_clim_start", "index_clim_end") := .(clim_start, clim_end)]

				indices_dt[.(plot, all_years[i]),
					c("year_start", "year_end") := .(all_years[i], all_years[i + 1])]
			}
		}
		print("100% done")
	}

	#### Create indices
	## For the state space model (*.stan file)
	indices = indices_state_space(trees_NFI = treeData)
	setkey(indices, plot_id, tree_id, year)

	## For the climate
	indices_climate(indices, climate, time_space)

	## For the average climate
	indices_avg = unique(indices[, .(plot_id, year, index_clim_start, index_clim_end)])
	setkey(indices_avg, plot_id, year)

	indices_climate_avg(indices_avg, climate)
	indices_avg = na.omit(indices_avg)

	## Add an index per plot
	indices[, plot_index := .GRP, by = plot_id]

	## Add an index per country
	indices[, nfi := stri_sub(plot_id, to = stri_locate_first(plot_id, regex = "_")[, "start"] - 1)]
	indices[, nfi_index := .GRP, by = nfi]

	## Add type (either parent or child, which correspond to primary and subsequent in the article)
	indices[, type := "child"]

	# Correct for those who are parent type
	indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]

	## Compute the number of growing years per individual
	indices[, nbYearsGrowth := max(year) - min(year), by = .(plot_id, tree_id)]

	## Compute the number of growth intervals per individual
	indices[, nbIntervalGrowth := .N - 1, by = .(plot_id, tree_id)] # -1 because between n points, there are only n - 1 intervals

	## Check-up
	checkUp = all(indices[, nbYearsGrowth <= index_clim_end - index_clim_start])
	if(!checkUp)
		stop("Suspicious indexing for nbYearsGrowth. Review the indices data.table")

	checkUp = all(indices_avg[, year_end - year_start == index_clim_end - index_clim_start])
	if(!checkUp)
		stop("Suspicious indexing, review the indices_avg data.table")

	return (list(indices = indices, indices_avg = indices_avg))
}

#### Predict for a given dataset
## Common variables
species = "Fagus sylvatica"
run = 1
tree_folder = paste0("./", species, "/")
if (!dir.exists(tree_folder))
	stop(paste0("Path not found for species <", species, ">."))

## Initial dbh
# Load tree data
tree_data = readRDS(paste0(tree_folder, "1_treeData.rds"))
# Take 5 individuals:
#	- A: 2 within a plot, measured twice only, with one negative growth
#	- B: 2 within an other plot, measured twice and three times, respectively
#	- C: 1 measured three times and that is between A and B (in the data)
# tree_data = tree_data[c(1:4, 23995:23997, 24641, 24642, 24664:24666)]
(n_indiv = unique(tree_data[, .(tree_id, plot_id)])[, .N])

# Get index initial dbh
parents_index = tree_data[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
last_child_index = tree_data[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]

if (length(parents_index) != n_indiv)
	stop("Dimension mismatch between parents_index and n_indiv")

dbh0 = tree_data[parents_index, dbh]
last_dbh = tree_data[last_child_index, dbh]

## Environmental data
# Get environmental time series
dataEnv_ls = getEnvSeries(species, run)
dataEnv_ls[["dataEnv"]][, row_id := seq_len(.N)]

## Indexing
# Space-time
treeCoords = unique(tree_data[, .(plot_id, x, y)])
index_min = unique(tree_data[tree_data[, .(m = .I[year == min(year)]), by = c("plot_id", "x", "y")]$m,
	.(plot_id, x, y, year)])
index_max = unique(tree_data[tree_data[, .(M = .I[year == max(year)]), by = c("plot_id", "x", "y")]$M,
	.(plot_id, x, y, year)])

treeCoords[index_min, on = "plot_id", min_year := i.year]
treeCoords[index_max, on = "plot_id", max_year := i.year]
treeCoords[, unique_year_id := paste(min_year, max_year, sep = "-")]

setkey(treeCoords, plot_id)

ts_length = tree_data[!is.na(growth), sum(deltaYear)] # Corresponds to the number of latent growth, named after tree rings data
if (ts_length != sum(tree_data[last_child_index, year] - tree_data[parents_index, year]))
	stop("Time series length is wrong!")

# Indices for environment (for start_ind and end_ind)
indices = indices_predict(treeData = tree_data, climate = dataEnv_ls[["dataEnv"]], time_space = treeCoords)

## Determine suspicious growth (either negative or 3*mean(growth))
meanGrowth = tree_data[!is.na(growth), mean(growth)]
suspicious_growth_index = which(tree_data[parents_index, growth] < 0 | tree_data[parents_index, growth] > 3*meanGrowth)
non_suspicious_growth_index = which(tree_data[parents_index, growth] >= 0 & tree_data[parents_index, growth] <= 3*meanGrowth)

#! WATCH OUT ON HOW THESE INDICES ARE USED, IT SHOULD BE SOMETHING LIKE:
#&	tree_data[parents_index[non_suspicious_growth_index]]
#! AND NOT DIRECTLY
#&	tree_data[non_suspicious_growth_index]

#! Check also deltaYear vs nbYearsGrowth! deltaYear is between two dbh measurements while nbYearsGrowth is between first and last measurements

## Predict
preds = predict(species = species, run = run, env = dataEnv_ls[["dataEnv"]], dbh0 = dbh0, suspicious = suspicious_growth_index,
	indices = indices[["indices"]], time_series_length = ts_length)

## Access results
generate_quantities_ssm = preds[["generate_quantities_ssm"]]
stanData_ssm = preds[["stanData_ssm"]]

generate_quantities_classic = preds[["generate_quantities_classic"]]
stanData_classic = preds[["stanData_classic"]]

# SSM -----------------
predicted_growth_ssm = stanData_ssm$sd_dbh * generate_quantities_ssm$draws("simulatedGrowth")
predicted_growth_avg_ssm = stanData_ssm$sd_dbh * generate_quantities_ssm$draws("simulatedGrowth_avg")
predicted_growth_avg_ssm_avg = apply(X = predicted_growth_avg_ssm, MARGIN = 3, FUN = mean)
predicted_growth_ssm_avg = apply(X = predicted_growth_ssm, MARGIN = 3, FUN = mean)

# hist(predicted_growth_ssm_avg - predicted_growth_avg_ssm_avg, prob = TRUE)
# (aa = mean(predicted_growth_ssm_avg - predicted_growth_avg_ssm_avg))
# (bb = sd(predicted_growth_ssm_avg - predicted_growth_avg_ssm_avg))
# curve(dnorm(x, aa, bb), add = TRUE, col = "#CD1A21", lwd = 3)

# predicted_dbh_ssm = stanData_ssm$sd_dbh * generate_quantities_ssm$draws("simulatedObservedDBH")
# predicted_dbh_avg_ssm = stanData_ssm$sd_dbh * generate_quantities_ssm$draws("simulatedObservedDBH_avg") # simulatedLatentDBH_avg
# predicted_dbh_ssm_avg = apply(X = predicted_dbh_ssm, MARGIN = 3, FUN = mean)
# predicted_dbh_avg_ssm_avg = apply(X = predicted_dbh_avg_ssm, MARGIN = 3, FUN = mean)

# Classic -------------
predicted_growth_classic = stanData_classic$sd_dbh * generate_quantities_classic$draws("simulatedGrowth")
predicted_growth_avg_classic = stanData_classic$sd_dbh * generate_quantities_classic$draws("simulatedGrowth_avg")
predicted_growth_avg_classic_avg = apply(X = predicted_growth_avg_classic, MARGIN = 3, FUN = mean)
predicted_growth_classic_avg = apply(X = predicted_growth_classic, MARGIN = 3, FUN = mean)

# predicted_dbh_classic = stanData_classic$sd_dbh * generate_quantities_classic$draws("simulatedObservedDBH")
# predicted_dbh_avg_classic = stanData_classic$sd_dbh * generate_quantities_classic$draws("simulatedObservedDBH_avg")
# predicted_dbh_classic_avg = apply(X = predicted_dbh_classic, MARGIN = 3, FUN = mean)
# predicted_dbh_avg_classic_avg = apply(X = predicted_dbh_avg_classic, MARGIN = 3, FUN = mean)


#### Reconstruct the average annual growth from the annual latent growth

# Reconstruct growth for simulations from latent growth
growth_dt = na.omit(tree_data)
latent_growth_ssm = predicted_growth_ssm
reconstructG = function(latent_growth_ssm, growth_dt)
{
	n_iter = dim(latent_growth_ssm)[1]
	n_chains = dim(latent_growth_ssm)[2]
	n_growth = growth_dt[, .N]

	residualsPred_ssm = posterior::as_draws_array(array(data = NA, dim = c(n_iter, n_chains, n_growth)))
	reconstructedG_ssm = posterior::as_draws_array(array(data = NA, dim = c(n_iter, n_chains, n_growth)))

	start = 1

	for (g in seq_len(n_growth))
	{
		nbYears = growth_dt[g, deltaYear]
		end = start + nbYears - 1

		reconstructedG_ssm[, , g] =
			posterior::as_draws_array(array(rowSums(latent_growth_ssm[, , start:end], dim = 2), dim = c(n_iter, n_chains, 1)))/nbYears

		residualsPred_ssm[, , g] =
			posterior::as_draws_array(array(rowSums(latent_growth_ssm[, , start:end], dim = 2), dim = c(n_iter, n_chains, 1)))/nbYears -
			growth_dt[g, growth]

		start = end + 1
		if (g %% 200 == 0)
			print(paste0(round(g*100/n_growth, 2), "% done"))
	}

	return(list(reconstructedG_ssm = reconstructedG_ssm, residualsPred_ssm = residualsPred_ssm))
}

avgAnnualGrowth = reconstructG(predicted_growth_ssm, growth_dt)

avgAnnualGrowth_avg = apply(X = avgAnnualGrowth[["reconstructedG_ssm"]], MARGIN = 3, FUN = mean)


avgAnnualGrowth_fit = residuals_fit(species, run)
reconstructedGrowth_ssm_avg = apply(X = avgAnnualGrowth_fit[["reconstructedGrowth_ssm"]], MARGIN = 3, FUN = mean)

plot(avgAnnualGrowth_avg ~ reconstructedGrowth_ssm_avg, pch = 19, col = "#44556644")

#### Plots
growth = growth_dt[, growth]

plot(avgAnnualGrowth_avg ~ growth, pch = 19, col = "#44556644")
plot(reconstructedGrowth_ssm_avg ~ growth, pch = 19, col = "#44556644")

ind = which(growth >= 0)
plot(reconstructedGrowth_ssm_avg[ind] ~ growth[ind], pch = 19, col = "#44556644")
(reg = lm(reconstructedGrowth_ssm_avg[ind] ~ 0 + growth[ind]))
abline(reg, col = "#CD1A21", lwd = 3)


lim_05 = quantile(growth, 0.05)
lim_95 = quantile(growth, 0.95)
ind = which(growth >= lim_05 & growth < lim_95)

plot(100*avgAnnualGrowth_fit[["residuals_ssm"]][ind]/growth[ind] ~ growth[ind], pch = 19, col = "#44556644")
abline(h = -20, col = "#CD1A21", lwd = 1, lty = 2)
abline(h = 0, col = "#CD1A21", lwd = 3)
abline(h = 20, col = "#CD1A21", lwd = 1, lty = 2)

plot(100*avgAnnualGrowth_fit[["residuals_classic"]][ind]/growth[ind] ~ growth[ind], pch = 19, col = "#44556644")
abline(h = -20, col = "#CD1A21", lwd = 1, lty = 2)
abline(h = 0, col = "#CD1A21", lwd = 3)
abline(h = 20, col = "#CD1A21", lwd = 1, lty = 2)



#### Compare the draws latent_growth with the simulated latent_growth
info_lastRun = getLastRun(path = tree_folder, begin = "^growth-", extension = "_main.rds$", format = "ymd", run = run, hour = TRUE)
ssm = readRDS(paste0(tree_folder, info_lastRun[["file"]]))
latent_G_ssm = ssm$draws("latent_growth")

dim(predicted_growth_ssm)
dim(latent_G_ssm)

aa = predicted_growth_ssm - latent_G_ssm
bb = apply(X = aa, MARGIN = 3, FUN = mean)

hist(bb, prob = TRUE)
mean(bb)
sd(bb)
curve(dlnorm(x, meanlog = log(mean(bb)) - 0.5*log(1 + var(bb)/mean(bb)^2), sdlog = sqrt(log(1 + var(bb)/(mean(bb))^2))),
	lwd = 2, col = "#CD1A21", add = TRUE)

curve(dlnorm(x, meanlog = log(mean(bb)) - 0.5*log(1 + var(bb)/mean(bb)^2), sdlog = sqrt(log(1 + var(bb)/(mean(bb))^2))),
	lwd = 2, col = "#CD1A21", to = 10, ylab = "logNormal")
abline(v = exp(log(mean(bb)) - 1.5*log(1 + var(bb)/mean(bb)^2)), col = "#332277", lwd = 2)
abline(v = exp(log(mean(bb)) + 0.5*log(1 + var(bb)/mean(bb)^2)), col = "#332277", lwd = 2, lty = 2)



latent_dbh = stanData_ssm$sd_dbh * generate_quantities_ssm$draws("simulatedLatentDBH_avg")
dim(latent_dbh)
zz = apply(X = latent_dbh[, , indices[["indices"]][type == "parent", index_gen]], MARGIN = 3, FUN = mean)

max(abs((zz - dbh0)/dbh0)) # Relative error, ok!

dbh0_latent = apply(X = stanData_ssm$sd_dbh*ssm$draws("latent_dbh_parents"), MARGIN = 3, FUN = mean)
hist(dbh0_latent - zz)

# ## Plots
# tree_data[last_child_index, index_by_country := seq_len(.N), by = nfi_id]
# index_lim = tree_data[, max(index_by_country, na.rm = TRUE), by = nfi_id][, V1]

# if (length(index_lim) > 1)
# {
# 	for (i in 2:length(index_lim))
# 		index_lim[i] = index_lim[i] + index_lim[i - 1]
# }

# ind = 1:max(index_lim)
# ind = ind[!(ind %in% suspicious_growth_index)]

# plot(predicted_dbh_avg_ssm_avg - tree_data[last_child_index, dbh], pch = 19)
# for (i in seq_along(index_lim))
# 	abline(v = index_lim[i], lwd = 2, col = "#CD212A")

# plot(predicted_dbh_avg_ssm_avg[ind] - tree_data[last_child_index][ind, dbh], pch = 19)

# plot(tree_data[last_child_index, dbh], predicted_dbh_avg_ssm_avg, pch = 19)
# points(tree_data[last_child_index][suspicious_growth_index, dbh], predicted_dbh_avg_ssm_avg[suspicious_growth_index], pch = 19, col = "#CD212A")

# plot(tree_data[last_child_index][ind, dbh], predicted_dbh_avg_ssm_avg[ind], pch = 19)

# mean(predicted_dbh_avg_ssm_avg - tree_data[last_child_index, dbh])
# plot(density(predicted_dbh_avg_ssm_avg - tree_data[last_child_index, dbh]))

# reg = lm(predicted_dbh_avg_ssm_avg ~ tree_data[last_child_index, dbh])
# reg = lm(predicted_dbh_avg_ssm_avg[ind] ~ tree_data[last_child_index][ind, dbh])


# resid = residuals_fit(species, run)

# reconstructedG = resid[["reconstructedGrowth_ssm"]]
# reconstructedG_avg = apply(X = reconstructedG, MARGIN = 3, FUN = mean)

# plot(reconstructedG_avg, tree_data[!is.na(growth), growth], pch = 19)

# plot(reconstructedG_avg[non_suspicious_growth_index] ~ tree_data[parents_index][non_suspicious_growth_index, growth], pch = 19,
# 	col = "#AABBCC88")
# reg = lm(reconstructedG_avg[non_suspicious_growth_index] ~ 0 + tree_data[parents_index][non_suspicious_growth_index, growth])
# abline(reg, col = "#CD1A21", lwd = 3)

# reg = lm(reconstructedG_avg ~ tree_data[!is.na(growth), growth])







# plot(reconstructedG_avg, reconstructedGPred_avg, pch = 19)
# reg = lm(reconstructedGPred_avg ~ 0 + reconstructedG_avg)
# abline(reg, col = "#CD1A21", lwd = 3)

# plot(reconstructedG_avg ~ tree_data[!is.na(growth), growth], pch = 19, col = "#22BBCC33", xlim = c(-15, 25))
# points(reconstructedGPred_avg ~ tree_data[!is.na(growth), growth], pch = 19, col = "#AABBCC45")
# points(sort(reconstructedGPred_avg) ~ sort(tree_data[!is.na(growth), growth]), pch = 19, col = "#212A7844")





# plot(predicted_dbh_avg_ssm_avg ~ tree_data[last_child_index, dbh], pch = 19)
# reg = lm(predicted_dbh_avg_ssm_avg ~ 0 + tree_data[last_child_index, dbh])
# abline(reg, col = "#CD1A21", lwd = 3)

# residuals = predicted_dbh_avg_ssm_avg - tree_data[last_child_index, dbh]
# plot(residuals ~ tree_data[last_child_index, dbh], pch = 19, col = "#01020122")
# reg = lm(residuals ~ tree_data[last_child_index, dbh])
# abline(reg, col = "#CD1A21", lwd = 3)




# pred_dbh_final_reconstructed = tree_data[!is.na(growth), dbh] + reconstructedG_avg
# plot(pred_dbh_final_reconstructed ~ tree_data[!is.na(growth), dbh], pch = 19)
# reg = lm(pred_dbh_final_reconstructed ~ 0 + tree_data[!is.na(growth), dbh])
# abline(reg, col = "#CD1A21", lwd = 3)

# residuals = rnorm(tree_data[!is.na(growth), .N], mean = pred_dbh_final_reconstructed, sd = 3) - tree_data[!is.na(growth), dbh]
# plot(residuals ~ tree_data[!is.na(growth), dbh], pch = 19, col = "#01020122")
# reg = lm(residuals ~ tree_data[!is.na(growth), dbh])
# abline(reg, col = "#CD1A21", lwd = 3)




# # In the case predicted_dbh_avg_ssm_avg is the latent dbh, not the observed, i.e., generate_quantities_ssm$draws("simulatedLatentDBH_avg")
# plot(predicted_dbh_avg_ssm_avg[indices[type == "child", index_gen]] ~ tree_data[!parents_index, dbh], pch = 19)
# reg = lm(predicted_dbh_avg_ssm_avg[indices[type == "child", index_gen]] ~ tree_data[!parents_index, dbh])
# abline(reg, col = "#CD1A21", lwd = 3)

# residuals = predicted_dbh_avg_ssm_avg[indices[type == "child", index_gen]] - tree_data[!parents_index, dbh]
# plot(residuals ~ tree_data[!is.na(growth), growth], pch = 19, col = "#01020122")
# reg = lm(residuals ~ tree_data[!is.na(growth), growth])
# abline(reg, col = "#CD1A21", lwd = 3)
