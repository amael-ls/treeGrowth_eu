
#### Aim of prog: To have a look on the French data and get some hints on the priors
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var
#

rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(cmdstanr)
library(DHARMa)
library(actuar) # For the Pareto distribution

options(max.print = 500)

#### Load data
trees = readRDS("../Inventories/FR IFN/trees_forest.rds")

trees[, dbh_5_in_mm := dbh_in_mm + increment_5_yrs]
trees[, ratio := dbh_5_in_mm/dbh_in_mm]

hist(trees[, ratio])
hist(trees[, dbh_in_mm])

if (!file.exists("./ratio.jpg"))
{
	jpeg("ratio.jpg", quality = 100, height = 1080, width = 1080)
	plot(trees[, dbh_in_mm], trees[, ratio], pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0",
		col = "#32323210")
	dev.off()
}

#### Quick stan model on the ratio
model = cmdstan_model("ratio_pareto_noMoment.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results_pareto = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results_pareto$cmdstan_diagnose()

# Alpha_0
plot_title = ggplot2::ggtitle("Posterior distribution alpha_0", "with medians and 80% intervals")
mcmc_areas(results_pareto$draws("alpha_0"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha_0")
mcmc_trace(results_pareto$draws("alpha_0")) + plot_title

# Alpha_1
plot_title = ggplot2::ggtitle("Posterior distribution alpha_1", "with medians and 80% intervals")
mcmc_areas(results_pareto$draws("alpha_1"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha_1")
mcmc_trace(results_pareto$draws("alpha_1")) + plot_title

#### Check-up correlations
cor(results_pareto$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0,
	results_pareto$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1)

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

alpha0_ptw = mean(results_pareto$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0)
alpha0_ptw_med = median(results_pareto$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0)

alpha1_ptw = mean(results_pareto$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1)
alpha1_ptw_med = median(results_pareto$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1)

fct_alpha = function(x, alpha0 = alpha0_ptw, alpha1 = alpha1_ptw)
	return (alpha0*x^alpha1)

avg = function(x, alpha0 = alpha0_ptw, alpha1 = alpha1_ptw)
	return (fct_alpha(x, alpha0, alpha1)/(fct_alpha(x, alpha0, alpha1) - 1));

n_rep = 250

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled_mean = numeric(n_rep * length(dbh)), sampled_med = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	alpha_i_mean = fct_alpha(dbh[i]);
	alpha_i_median = fct_alpha(dbh[i], alpha0_ptw_med, alpha1_ptw_med);
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled_mean := rpareto1(n_rep, shape = alpha_i_mean, min = 1)]
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled_med := rpareto1(n_rep, shape = alpha_i_median, min = 1)]
}

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2, add = TRUE)
if (n_rep < 300)
	points(dt[, rep_dbh], dt[, sampled_mean], pch = 19, col = "#4DA6FF02")

#### Check variance per dbh class
diam_class = seq(min(dbh) - 1, max(dbh) + 1, by = 10)
var_class = numeric(length(diam_class))
var_class_sim = numeric(length(diam_class))
for (i in 1:(length(diam_class) - 1))
{
	var_class[i] = var(trees[dbh_in_mm >= diam_class[i] & dbh_in_mm < diam_class[i + 1], ratio])
	var_class_sim[i] = var(dt[rep_dbh >= diam_class[i] & rep_dbh < diam_class[i + 1], sampled_mean])
}

plot(diam_class, var_class - var_class_sim, pch = 19, col = "#CD212A", ylab = "difference data - simulations")

#### Check residuals (using DHARMa)
## Reshape simulations into a matrix length(dbh) x n_rep
sims_mean = matrix(data = dt[, sampled_mean], nrow = n_rep, ncol = length(dbh)) # each column is for one data point (the way I organised the dt)
sims_mean = t(sims_mean) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims_mean, observedResponse = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio]) #,
	# fittedPredictedResponse = avg(dbh, alpha0_ptw_med, alpha1_ptw_med))
plot(forDharma)

#?#################################################################
#*###########    MODEL USING EITZEL 2013'S FORMULA    #############
#?#################################################################
#### Quick stan model on the ratio
model = cmdstan_model("eitzel_dbh.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	dbh0 = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm],
	dbh1 = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_5_in_mm])

results_eitzel = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, # init = init_lambda,
	max_treedepth = 11, adapt_delta = 0.95)

results_eitzel$cmdstan_diagnose()

# Alpha
plot_title = ggplot2::ggtitle("Posterior distribution alpha", "with medians and 80% intervals")
mcmc_areas(results_eitzel$draws("alpha"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha")
mcmc_trace(results_eitzel$draws("alpha")) + plot_title

# Beta
plot_title = ggplot2::ggtitle("Posterior distribution beta", "with medians and 80% intervals")
mcmc_areas(results_eitzel$draws("beta"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for beta")
mcmc_trace(results_eitzel$draws("beta")) + plot_title

# Sigma_eps
plot_title = ggplot2::ggtitle("Posterior distribution sigma_eps", "with medians and 80% intervals")
mcmc_areas(results_eitzel$draws("sigma_eps"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for sigma_eps")
mcmc_trace(results_eitzel$draws("sigma_eps")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

alpha_ptw = mean(results_eitzel$draws(variables = c("alpha"), format = "draws_df")$alpha)
alpha_ptw_med = median(results_eitzel$draws(variables = c("alpha"), format = "draws_df")$alpha)

beta_ptw = mean(results_eitzel$draws(variables = c("beta"), format = "draws_df")$beta)
beta_ptw_med = median(results_eitzel$draws(variables = c("beta"), format = "draws_df")$beta)

sigma_eps_ptw = mean(results_eitzel$draws(variables = c("sigma_eps"), format = "draws_df")$sigma_eps)
sigma_eps_ptw_med = median(results_eitzel$draws(variables = c("sigma_eps"), format = "draws_df")$sigma_eps)

avg = function(x, alpha = alpha_ptw, beta = beta_ptw)
	return (alpha + beta*x);

n_rep = 250

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	avg_i = avg(dbh[i]);
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rnorm(n = n_rep, mean = avg_i, sd = sigma_eps_ptw)]
}

dt[, ratio := sampled/rep_dbh]

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_5_in_mm], # ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2, add = TRUE)
if (n_rep < 300)
	points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF01")

diam_class = seq(min(dbh) - 1, max(dbh) + 1, by = 10)
var_class = numeric(length(diam_class))
var_class_sim = numeric(length(diam_class))
for (i in 1:(length(diam_class) - 1))
{
	var_class[i] = var(trees[dbh_in_mm >= diam_class[i] & dbh_in_mm < diam_class[i + 1], ratio])
	var_class_sim[i] = var(dt[rep_dbh >= diam_class[i] & rep_dbh < diam_class[i + 1], sampled])
}

plot(diam_class, var_class - var_class_sim, pch = 19, col = "#CD212A")

#### Check residuals (using DHARMa)
## Reshape simulations into a matrix length(dbh) x n_rep
sims = matrix(data = dt[, sampled], nrow = n_rep, ncol = length(dbh)) # each column is for one data point (the way I organised the dt)
sims = t(sims) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims,
	observedResponse = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_5_in_mm],
	fittedPredictedResponse = avg(dbh, alpha0_ptw_med, alpha1_ptw_med))
plot(forDharma)
