
#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#### Tool functions
## Function to compute all the combinations of overlaps for one parameter
overlap_fct = function(polygons_ls, n_runs)
{
	checkPolygons = sapply(polygons_ls, terra::is.valid)

	if (!all(checkPolygons))
	{
		nn = names(polygons_ls)
		warning(paste("The polygons", nn[!checkPolygons], "were not valid, but have been forced to be via terra::makeValid"))
		for (polyg in nn[!checkPolygons])
			polygons_ls[[polyg]] = terra::makeValid(polygons_ls[[polyg]])
	}

	overlap_ls = vector(mode = "list", n_runs - 1)
	names(overlap_ls) = paste0(2:n_runs, "_elements")

	for (nbElements in 2:n_runs) # Number of elements in the combination
	{
		combinations_ls = combn(x = 1:n_runs, m = nbElements, simplify = TRUE)
		results_dt = data.table(combination = character(ncol(combinations_ls)), overlap = numeric(ncol(combinations_ls)))
		
		for (j in seq_len(results_dt[, .N])) # Loop among the combinations with 'nbElements' runs
		{
			selected_runs = paste0("run_", combinations_ls[, j])
			currentCombination = paste(combinations_ls[, j], collapse = "_")

			union_polygon = polygons_ls[[selected_runs[1]]]
			inter_polygon = polygons_ls[[selected_runs[1]]]

			for (currentRun in selected_runs[2:length(selected_runs)]) # Loop among the selected runs
			{
				union_polygon = terra::union(union_polygon, polygons_ls[[currentRun]])
				union_polygon = terra::aggregate(union_polygon)
				inter_polygon = terra::intersect(inter_polygon, polygons_ls[[currentRun]])
			}
			results_dt[j, combination := currentCombination]
			area = suppressWarnings(ifelse(is.polygons(inter_polygon), terra::expanse(inter_polygon), 0))
			suppressWarnings(results_dt[j, overlap := area/terra::expanse(union_polygon)])
		}
		overlap_ls[[paste0(nbElements, "_elements")]] = results_dt
	}

	overlap_dt = rbindlist(overlap_ls, idcol = "elements")
	overlap_dt[, elements := stri_sub(str = elements, to = stri_locate_first(str = elements, regex = "_")[, "start"] - 1)]
	
	return(overlap_dt)
}

## Function to plot growth posteriors for all runs and compute their overlaps
lazyPosteriorGrowth = function(draws, printPlot = FALSE)
{
	n_runs = length(draws)

	# Get posterior
	density_from_draws = vector(mode = "list", length = n_runs)
	x = vector(mode = "list", length = n_runs)
	y = vector(mode = "list", length = n_runs)
	polygons_ls = vector(mode = "list", length = n_runs)
	names(polygons_ls) = paste0("run_", 1:n_runs)

	for (run in seq_len(n_runs))
	{
		density_from_draws[[run]] = density(sd_dbh*draws[[run]], n = 512)
		x[[run]] = density_from_draws[[run]]$x
		y[[run]] = density_from_draws[[run]]$y

		temporary = cbind(id = 1, part = 1, x[[run]], y[[run]])
		polygons_ls[[run]] = vect(temporary, type = "polygons")
	}
	min_x = min(sapply(x, min))
	max_x = max(sapply(x, max))
	max_y = max(sapply(y, max))

	min_x = ifelse(min_x < 0, 1.1*min_x, 0.9*min_x) # To extend 10% from min_x
	max_x = ifelse(max_x < 0, 0.9*max_x, 1.1*max_x) # To extend 10% from max_x
	
	# Plot posterior and compute overlap
	if (printPlot)
	{
		colours = MetBrewer::met.brewer("Hokusai3", n_runs)
		colours_str = grDevices::colorRampPalette(colours)(n_runs)
		colours_str_pol = paste0(colours_str, "66")
		plot(0, type = "n", xlim = c(min_x, max_x), ylim = c(0, max_y), ylab = "frequence", main = "", xlab = "")
		for (i in seq_len(n_runs))
		{
			lines(x = density_from_draws[[i]]$x, y = density_from_draws[[i]]$y, col = colours_str[i], lwd = 2)
			polygon(density_from_draws[[i]], col = colours_str_pol[i])
		}
	}

	overlap = overlap_fct(polygons_ls = polygons_ls, n_runs = n_runs)[elements == n_runs, overlap]
	return (overlap)
}

## Other functions
source("./toolFunctions.R")

#### Common variables
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
	stop("Supply the species", call. = FALSE)

species = as.character(args[1])
path = paste0("./", species, "/")

print(paste("Running for species", species))

n_runs = 4 # Number of runs used in growth_subsample.R

## Species informations
infoSpecies = readRDS("./speciesInformations.rds")

threshold_indiv = 12000 # Minimal number of individuals required to use multi runs
threshold_time = as.Date("2023/01/01") # Results anterior to this date will not be considered

infoSpecies[, multiRun := if (n_indiv > threshold_indiv) TRUE else FALSE, by = speciesName_sci]
infoSpecies[, processed := isProcessed(path = speciesName_sci, multi = multiRun, lim_time =  threshold_time,
	extension = "_de-fr-sw_12000_main.rds$", lower = 1, upper = n_runs), by = speciesName_sci]
infoSpecies = infoSpecies[(processed)]

if (!(species %in% infoSpecies[, speciesName_sci]))
	stop(paste0("Species <", species, "> not processed"))

if (infoSpecies[species, n_indiv] < threshold_indiv)
	stop("No need to compute overlap, there is only one run")

nb_nfi = infoSpecies[species, n_nfi]

## Parameters
params = c("averageGrowth", "dbh_slope", "dbh_slope2", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
	"ph_slope", "ph_slope2", "competition_slope", "sigmaProc", "etaObs")

if (nb_nfi > 1)
	params = expand(params, nb_nfi)[["new_names"]]

n_params = length(params)

## Posteriors object
polygons_ls = vector(mode = "list", length = n_params)
names(polygons_ls) = params

for (currentParam in params)
{
	polygons_ls[[currentParam]] = vector(mode = "list", length = n_runs)
	names(polygons_ls[[currentParam]]) = paste0("run_", 1:n_runs)
}

#### Compute overlap posteriors
## Load results and extract posteriors
print("Starting loading results")
for (i in 1:n_runs)
{
	info_lastRun = getLastRun(path = path, extension = "_main.rds$", run = i)
	lastRun = info_lastRun[["file"]]
	results = readRDS(paste0(path, lastRun))

	currentRun = paste0("run_", i)

	for (currentParam in params)
	{
		currentDensity = density(results$draws(currentParam), n = 2048)
		x = currentDensity$x
		y = currentDensity$y

		temporary = cbind(id = 1, part = 1, x, y)
		polygons_ls[[currentParam]][[currentRun]] = vect(temporary, type = "polygons")
	}
	
	print(paste("Run", i, "done"))
}

## Overlaps for all the combinations
overlap_ls = vector(mode = "list", length(params))
names(overlap_ls) = params

for (currentParam in params)
	overlap_ls[[currentParam]] = overlap_fct(polygons_ls[[currentParam]], n_runs)

overlap_dt = rbindlist(overlap_ls, idcol = "parameter")

saveRDS(overlap_dt, paste0(path, "overlap_dt.rds"))

print(paste(species, "done"))


####! CRASH TEST ZONE -------------------------------------------------------
set.seed(123)
n = 20000

e1 = rnorm(n, mean = 0, sd = 10)
e2 = runif(n, min = -3, max = 100)
dbh = runif(n, min = 17, max = 1000)

beta0 = 0.5
beta1 = 3
beta2 = -1

gamma1 = -1.5
gamma2 = 0.78

delta1 = -0.7
delta2 = -0.08

sd_dbh = sd(dbh)

m_e1 = mean(e1)
m_e2 = mean(e2)

s_e1 = sd(e1)
s_e2 = sd(e2)

e1_tilde = (e1 - m_e1)/s_e1
e2_tilde = (e2 - m_e2)/s_e2
dbh_tilde = dbh/sd_dbh

meanlog = beta0 + beta1*dbh_tilde + beta2*dbh_tilde^2 + gamma1*e1_tilde + gamma2*e2_tilde + delta1*e1_tilde^2 + delta2*e2_tilde^2
sdlog = 0.15

G_tilde = rlnorm(n, meanlog = meanlog, sdlog = sdlog)
G = sd_dbh * G_tilde

meanlog_rescale = beta0 + log(sd_dbh) - (m_e1/s_e1*gamma1 + m_e2/s_e2*gamma2) + (m_e1^2/s_e1^2*delta1 + m_e2^2/s_e2^2*delta2) +
	beta1/sd_dbh * dbh + beta2/sd_dbh^2 * dbh^2 +
	(gamma1/s_e1 - 2*delta1*m_e1/s_e1^2)*e1 + (gamma2/s_e2 - 2*delta2*m_e2/s_e2^2)*e2 +
	delta1/s_e1^2*e1^2 + delta2/s_e2^2*e2^2

pred_tilde = exp(meanlog + sdlog^2/2)
pred = exp(meanlog_rescale + sdlog^2/2)
pred2 = sd_dbh*exp(meanlog + sdlog^2/2)

all.equal(pred2, pred)

####! END CRASH TEST ZONE ---------------------------------------------------
