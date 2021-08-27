
# This prog generates data for toy.stan, toyPara_GPUs.stan and toyPara_reduce_sum.stan
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var

####
rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(cmdstanr)
library(shiny)

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
record_min = 3
record_max = 6

n_census_per_plot = sample(x = record_min:record_max, size = n_plots, replace = TRUE)

if (any(n_census_per_plot < 2))
	stop("There should be at least 2 census per plot")

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
intercept = 0 # 2.4

slope_precip = 0 # 0.004
quad_slope_precip = 0 # -0.00001

slope_dbh = 1.1

sigma_process = 10
sigma_measure = 4

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
	initial_dbh = rgamma(n_indiv_per_plot[xy], shape = 150^2/500, rate = 150/500)
	index = seq(from = 1, to = nb_measures_plot, by = n_census_per_plot[xy]) + count
	# treeStates_dt[index, true_dbh := initial_dbh]
	treeStates_dt[index, true_dbh := rnorm(n = n_indiv_per_plot[xy], mean = initial_dbh, sd = sigma_process)]

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
setorder(treeStates_dt, plot_id, tree_id, year)
treeStates_dt[, unique_id := 1:.N]

#### Generate partial data set (keep only first and last measurements)
kept_rows = sort(c(treeStates_dt[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1],
	treeStates_dt[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]))

##! Quick and dirty solution to also keep an intermediate measure between year_min and year_max
cheating = integer(sum(n_indiv_per_plot))
for (i in 1:sum(n_indiv_per_plot))
	cheating[i] = round((kept_rows[2*i - 1] + kept_rows[2*i])/2)


if (any(cheating %in% kept_rows))
	warning("Some individuals did not get an extra measurment")

kept_rows = sort(c(kept_rows, cheating))

##! Quick and dirty solution to keep everything
kept_rows = unique(1:treeStates_dt[, .N])

treeData = treeStates_dt[kept_rows]
print(paste0(round(treeData[, .N]*100/treeStates_dt[, .N], 2), " % of the data kept"))

setorder(treeData, plot_id, tree_id, year)

#### plots
## General plot
plot(treeStates_dt[, unique_id], treeStates_dt[, true_dbh], pch = 19, cex = 0.7, col = "red", ty = "o", xlab = "unique id", ylab = expression(observed_dbh[i](t)),
	ylim = c(min(treeStates_dt[, .(true_dbh, observed_dbh)]) - abs(min(treeStates_dt[, .(true_dbh, observed_dbh)]))/10,
		max(treeStates_dt[, .(true_dbh, observed_dbh)] + max(treeStates_dt[, .(true_dbh, observed_dbh)]/10))), las = 1)
points(treeData[, unique_id], treeData[, observed_dbh], pch = 3, cex = 0.8, col = "blue", ty = "o", lty = 3)

legend("top", legend = c("Obs. (with missing data)", "True states"), pch = c(3, 19),
	col = c("blue", "red"), lty = c(3, 1), horiz=TRUE, bty="n", cex=0.9)

## Few plots
par(mfrow = c(3, 1))
for (i in 1:3)
{
	plot(treeStates_dt[(tree_id == i) & (plot_id == 1), unique_id], treeStates_dt[(tree_id == i) & (plot_id == 1), true_dbh],
		pch = 19, cex = 0.7, col = "red", ty = "o", xlab = "unique id", ylab = expression(observed_dbh[i](t)),
		las = 1)
	points(treeData[(tree_id == i) & (plot_id == 1), unique_id],
		treeData[(tree_id == i) & (plot_id == 1), observed_dbh],
		pch = 3, cex = 0.8, col = "blue", ty = "o", lty = 3)

	legend("top", legend = c("Obs. (with missing data)", "True states"), pch = c(3, 19),
		col = c("blue", "red"), lty = c(3, 1), horiz=TRUE, bty="n", cex=0.9)
}
dev.off()

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



###########################################################?
######## 		THIRD PART: Run stan program 		########
###########################################################?

#### Tool function
## Initiate Y_gen with reasonable value (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh_parents", "years_indiv", "average_G", "n_hiddenState")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide chain_id, Yobs, and years_indiv")
	
	# print(providedArgs)

	# chain_id = providedArgs[["chain_id"]]
	dbh_parents = providedArgs[["dbh_parents"]]
	years_indiv = providedArgs[["years_indiv"]]
	average_G = providedArgs[["average_G"]]
	n_hiddenState = providedArgs[["n_hiddenState"]]

	Y_gen = numeric(n_hiddenState)

	count = 0

	for (i in 1:n_indiv) # Not that this forbid trees to shrink
	{
		Y_gen[count + 1] = rgamma(1, shape = dbh_parents[i]^2, rate = dbh_parents[i]) # Mean = dbh_parents[i], Variance = 1
		for (j in 2:years_indiv[i])
			Y_gen[count + j] = Y_gen[count + j - 1] + average_G[i] + rgamma(1, shape = 0.25, rate = 0.5)

		count = count + years_indiv[i];
	}

	return(list(Y_generated = Y_gen))
}

## A hint on how the priors look like
# Define UI for sliders
ui = fluidPage(
	# ---- App title
	titlePanel("Priors preview"),

	# ---- Sidebar layout with input and output definitions
	sidebarLayout(
		# ---- Sidebar to demonstrate various slider options
		sidebarPanel(
			# ---- Input: Mean associated to normal distribution plot
			sliderInput("mean_dnorm", "mean (normal distrib):",
				min = -50, max = 50,
				value = 0, step = 0.1),

			# ---- Input: Std deviation associated to normal distribution plot
			sliderInput("sd_dnorm", "standard deviation (normal distrib):",
				min = 0, max = 20,
				value = 10, step = 0.01),

			# ---- Input: Range to plot normal distribution
			sliderInput("range_dnorm", "Range normal:",
				min = -80, max = 80,
				value = c(-40, 40)),
			
			# ---- Input: Mean associated to normal distribution plot
			sliderInput("mean_dgamma", "mean (gamma distrib):",
				min = 0, max = 200,
				value = 5, step = 1),

			# ---- Input: Std deviation associated to normal distribution plot
			sliderInput("var_dgamma", "variance (gamma distrib):",
				min = 0, max = 1000,
				value = 5, step = 1),

			# ---- Input: Range to plot gamma distribution
			sliderInput("range_dgamma", "Range gamma:",
				min = 0, max = 200,
				value = c(0, 20)),

			# ---- Input: Integration range
			sliderInput("integration_range", "a: P(X < a):",
				min = 0, max = 200,
				value = 1, step = 0.5)
		),

		# ---- Main panel for displaying outputs
		mainPanel(
			# ---- Output: Plot of the distributions and display shape/scale for gamma distrib
			tableOutput("values"),
			plotOutput(outputId = "normal_pdf"),
			plotOutput(outputId = "gamma_pdf")
		)
	)
)

## Define server logic for slider examples
server = function(input, output) {
	# --- Reactive expression to display the value of shape and rate for gamma distrib
	sliderValues = reactive({
	data.frame(
	Name = c("Gamma distribution, shape:", "Gamma distribution, rate:", "P(x < a), gamma distrib:"),
	Value = as.character(c(
		paste(input$mean_dgamma, "^2/", input$var_dgamma, sep = ""),
		paste(input$mean_dgamma, "/", input$var_dgamma, sep = ""),
		integrate(
			function(x){return (dgamma(x, shape = input$mean_dgamma^2/input$var_dgamma, rate = input$mean_dgamma/input$var_dgamma))},
			lower = 0, upper = input$integration_range)$value)),
	stringsAsFactors = FALSE)
	})

	# ---- Show the values in an HTML table
	output$values <- renderTable({sliderValues()})

	output$normal_pdf = renderPlot({curve(dnorm(x, input$mean_dnorm, input$sd_dnorm),
		input$range_dnorm[1], input$range_dnorm[2], lwd = 5, col = "#0072B5",
		ylab = "Normal pdf")})

	output$gamma_pdf = renderPlot({
		curve(dgamma(x,
			shape = input$mean_dgamma^2/input$var_dgamma, rate = input$mean_dgamma/input$var_dgamma),
			input$range_dgamma[1], input$range_dgamma[2], lwd = 5, col = "#363945",
			ylab = "Gamma pdf")

		sh <<- input$mean_dgamma^2/input$var_dgamma # https://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
		ra <<- input$mean_dgamma/input$var_dgamma

		DescTools::Shade(
			dgamma(x, shape = sh, rate = ra),
			breaks = c(0.01, input$integration_range), col = "#A1A1A166", density = NA)
		})
}

# shinyApp(ui, server)

#### Load data
## Data are already in memory
precip = precipitations_yearly
rm(precipitations_yearly)

# Set everyone to child type
indices[, type := "child"]

# Correct for those who are parent type
indices[indices[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1], type := "parent"]

# Compute nb years per individual
indices[, nbYearsPerIndiv := max(year) - min(year) + 1, by = .(plot_id, tree_id)]

checkUp = all(indices[, nbYearsPerIndiv == index_precip_end - index_precip_start + 1])
if(!checkUp)
	stop("Suspicious indexing. Review the indices data.table")

if (length(precip) != indices[.N, index_precip_end])
	stop("Dimensions mismatch between climate and indices")

if (indices[, .N] != treeData[, .N])
	stop("Dimensions mismatch between indices and treeData")

n_obs = treeData[, .N]
print(paste("Number of data:", n_obs))

n_indiv = unique(treeData[, .(tree_id, plot_id)])[, .N]
print(paste("Number of individuals:", n_indiv))

nbYearsPerIndiv = unique(indices[, .(tree_id, plot_id, nbYearsPerIndiv)])[, nbYearsPerIndiv]
if (length(nbYearsPerIndiv) != n_indiv)
	stop("Dimension mismatch between nbYearsPerIndiv and n_indiv")

parents_index = treeData[, .I[which.min(year)], by = .(plot_id, tree_id)][, V1]
children_index = treeData[, .I[which(year != min(year))], by = .(plot_id, tree_id)][, V1]
last_child_index = treeData[, .I[which.max(year)], by = .(plot_id, tree_id)][, V1]
not_parent_index = 1:indices[.N, index_gen]
not_parent_index = not_parent_index[!(not_parent_index %in% indices[type == "parent", index_gen])]

if (length(parents_index) != n_indiv)
	stop("Dimensions mismatch between parents_index and n_indiv")

if (length(children_index) != n_obs - n_indiv)
	stop("Dimensions mismatch between children_index and number of children")

if (length(not_parent_index) != indices[.N, index_gen] - n_indiv)
	stop("Dimensions mismatch between not_parent_index, n_hiddenState, and n_indiv")

Y_generated_0 = rnorm(n_indiv, treeData[parents_index, observed_dbh], 5)

#### Stan model
## Define stan variables
# Common variables
maxIter = 6e3
n_chains = 3

# Data to provide
stanData = list(
	# Number of data
	n_indiv = n_indiv, # Number of individuals
	n_precip = length(precip), # Dimension of the climate vector
	n_obs = n_obs, # Number of trees observations
	n_hiddenState = indices[.N, index_gen], # Dimension of the state space vector
	n_children = n_obs - n_indiv, # Number of children trees observations = n_obs - n_indiv
	nbYearsPerIndiv = nbYearsPerIndiv, # Number of years for each individual

	# Indices
	parents_index = parents_index, # Index of each parent in the observed data
	parentsObs_index = indices[type == "parent", index_gen], # Corresponding index of observed parents in Y_generated
	children_index = children_index, # Index of children in the observed data
	childrenObs_index = indices[type == "child", index_gen], # Corresponding index of observed children in Y_generated
	climate_index = indices[type == "parent", index_precip_start], # Index of the climate associated to each parent
	not_parent_index = not_parent_index, # Index in Y_generated of states that cannot be compared to data

	# Initial state, used to start the states and compare with Y_generated[parentsObs_index] (i.e., a prior parameter)
	Y_generated_0 = Y_generated_0,

	# Observations
	Yobs = treeData[, observed_dbh],

	# Explanatory variable
	precip = precip, # Precipitations

	# Parameter for parralel calculus
	grainsize = 1
)

# Initial value for states only
average_G = (treeData[last_child_index, observed_dbh] - treeData[parents_index, observed_dbh])/
	(treeData[last_child_index, year] - treeData[parents_index, year])

if (length(average_G) != n_indiv)
	stop("Dimensions mismatch between average_G and number of individuals")

initVal_Y_gen = lapply(1:n_chains, init_fun, dbh_parents = treeData[parents_index, observed_dbh],
 	years_indiv = nbYearsPerIndiv, average_G = average_G,
 	n_hiddenState = indices[.N, index_gen])

length(initVal_Y_gen)

for (i in 1:n_chains)
	print(range(initVal_Y_gen[[i]]))

## Compile stan model
# model = stan_model(file = "./toy.stan")
# model = cmdstan_model("toyPara_reduce_sum.stan", cpp_options = list(stan_threads = TRUE)) # list(stan_threads = TRUE, stan_opencl = TRUE)
# model = cmdstan_model("toyPara_GPUs.stan", cpp_options = list(stan_opencl = TRUE))

model = cmdstan_model("toyPara_GPUs.stan")

## Run model
start = proc.time()

# results = stan(file = "toy.stan", data = stanData, cores = n_chains,
# 	iter = maxIter, chains = n_chains, init = initVal_Y_gen)

# results = model$sample(data = stanData, parallel_chains = n_chains, threads_per_chain = 8, refresh = 2,
# 	iter_warmup = maxIter/2, iter_sampling = maxIter/2, chains = n_chains, init = initVal_Y_gen)

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains, # threads_per_chain = 2,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, init = initVal_Y_gen, max_treedepth = 14, adapt_delta = 0.95)

# results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 2, chains = n_chains,
# 	opencl_ids = c(0, 0), iter_warmup = maxIter/2, iter_sampling = maxIter/2, init = initVal_Y_gen)

proc.time() - start

# results$save_object(file = paste0("./toyPara_GPUs_nonInformative_noProcessError_", format(Sys.time(), "%d-%m-%Y_%Hh%M"), ".rds"))

results$cmdstan_diagnose()

#### Plots of posterior distributions and traces
# ## Intercept
# plot_title = ggplot2::ggtitle("Traces for intercept")
# mcmc_trace(results$draws("intercepts")) + plot_title

# plot_title = ggplot2::ggtitle("Posterior distribution intercept", "with medians and 80% intervals")
# mcmc_areas(results$draws("intercepts"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = intercept, color = "#FFA500")

## Slope for dbh
plot_title = ggplot2::ggtitle("Posterior distribution slope dbh", "with medians and 80% intervals")
mcmc_areas(results$draws("slopes_dbh"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = slope_dbh, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for slope dbh")
mcmc_trace(results$draws("slopes_dbh")) + plot_title

# ## Slope for precipitation
# plot_title = ggplot2::ggtitle("Posterior distribution slope precip", "with medians and 80% intervals")
# mcmc_areas(results$draws("slopes_precip"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = slope_precip, color = "#FFA500")

# plot_title = ggplot2::ggtitle("Traces for slope precip")
# mcmc_trace(results$draws("slopes_precip")) + plot_title

# ## Slope for precipitation (quadratique term)
# plot_title = ggplot2::ggtitle("Posterior distribution slope precip (quadratic term)", "with medians and 80% intervals")
# mcmc_areas(results$draws("quad_slopes_precip"), prob = 0.8) + plot_title +
# 	ggplot2::geom_vline(xintercept = quad_slope_precip, color = "#FFA500")

# plot_title = ggplot2::ggtitle("Traces for slope precip (quadratic term)")
# mcmc_trace(results$draws("quad_slopes_precip")) + plot_title

## Measurement error
plot_title = ggplot2::ggtitle("Posterior distribution measure error", "with medians and 80% intervals")
mcmc_areas(results$draws("measureError"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigma_measure, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for measure error")
mcmc_trace(results$draws("measureError")) + plot_title

## Process error
plot_title = ggplot2::ggtitle("Traces for process error")
mcmc_trace(results$draws("processError")) + plot_title

plot_title = ggplot2::ggtitle("Posterior distribution process error", "with medians and 80% intervals")
mcmc_areas(results$draws("processError"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigma_process, color = "#FFA500")


# intercept = 2.4

# slope_precip = 0.004
# quad_slope_precip = -0.00001

# slope_dbh = 1.1

# sigma_process = 4
# sigma_measure = 1

## Check-up states
chosen_state = 50
if (chosen_state > treeStates_dt[, .N])
	stop(paste("chosen_state must be smaller than", treeStates_dt[, .N]))

plot_title = ggplot2::ggtitle(paste("Traces for state", chosen_state))
mcmc_trace(results$draws(paste0("Y_generated[", chosen_state, "]"))) + plot_title

plot_title = ggplot2::ggtitle(paste("Posterior distribution for state", chosen_state), "with medians and 80% intervals")
mcmc_areas(results$draws(paste0("Y_generated[", chosen_state, "]")), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = treeStates_dt[chosen_state, true_dbh], color = "#FFA500")

# Colours from https://www.w3schools.com/colors/colors_2021.asp

# classic_blue = "#34568B"
# flame_scarlet = "#CD212A"
# saffron = "#FFA500"
# biscay_green = "#56C6A9"
# chive = "#4B5335"
# faded_denim = "#798EA4"
# orange_peel = "#FA7A35"
# mosaic_blue = "#00758F"
# sunlight = "#EDD59E"
# coral_pink = "#E8A798"
# cinnamon_stick = "#9C4722"
# grape_compote = "#6B5876"
# lark = "#B89B72"
# navy_blazer = "#282D3C"
# brilliant_white = "#EDF1FF"
# ash = "#A09998"
# amberglow = "#DC793E"
# samba = "#A2242F"
# sandstone = "#C48A69"
# green_sheen = "#D9CE52"
# rose_tan = "#D19C97"
# ultramarine_green = "#006B54"
# fired_brick = "#6A2E2A"
# peach_nougat = "#E6AF91"
# magenta_purple = "#6C244C"

# marigold = "#FDAC53"
# cerulean = "#9BB7D4"
# rust = "#B55A30"
# illuminating = "#F5DF4D"
# french_blue = "#0072B5"
# green_ash = "#A0DAA9"
# burnt_coral = "#E9897E"
# mint = "#00A170"
# amethyst_orchid = "#926AA6"
# raspberry_sorbet = "#D2386C"
# inkwell = "#363945"
# ultimate_gray = "#939597"
# buttercream = "#EFE1CE"
# desert_mist = "#E0B589"
# willow = "#9A8B4F"

