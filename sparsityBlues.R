
## This script is to follow the blog post from Michael Bettancourt, that can be found at:
# https://betanalpha.github.io/assets/case_studies/modeling_sparsity.html#1_Fading_into_Irrelevance

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
source("./lazyTools.R")

#### Generate data
## Common variables
c_light = c("#DCBCBC")
c_light_highlight = c("#C79999")
c_mid = c("#B97C7C")
c_mid_highlight = c("#A25050")
c_dark = c("#8F2727")
c_dark_highlight = c("#7C0000")

c_light_trans = c("#DCBCBC80")
c_dark_trans = c("#8F272780")
c_green_trans = c("#00FF0080")

set.seed(58533858)

## data
K = 9
theta_true = c(-0.09604655,	0.02063505, 0.01246716, -0.04487767, -0.13519031, 0.09421304, -29.449233, 32.997172, 18.517443)

N = K
sigma = 0.5
context_idx = rep(1:K, 1)

y = rnorm(K, theta_true[context_idx], sigma)

data = list("N" = N, "K" = K, "context_idx" = context_idx, "y" = y, "sigma" = sigma)

#### Section 2.1: The Limitations of Normal Population Models
## Common data
n_chains = 4

## Fit data (narrow)
model_sec2.1_narrow = cmdstan_model("./normal_narrow.stan")

results = model_sec2.1_narrow$sample(data = data, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = 1000, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 13, adapt_delta = 0.8)

results$cmdstan_diagnose()

samples = results$draws("theta")

par(mfrow=c(3, 3))

for (k in 1:9)
{
	hist(samples[, , k], breaks = seq(-35, 35, 0.25), main = paste("k = ", k), col = c_dark, border = c_dark_highlight,
		xlab = "theta", yaxt = 'n', ylab = "")
	abline(v = theta_true[k], col = "white", lwd = 2)
	abline(v = theta_true[k], col = "black", lwd = 1)
}

## Fit data (wide)
model_sec2.1_wide = cmdstan_model("./normal_wide.stan")

results = model_sec2.1_wide$sample(data = data, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = 1000, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 13, adapt_delta = 0.8)

results$cmdstan_diagnose()

samples = results$draws("theta")

par(mfrow=c(3, 3))

for (k in 1:9)
{
	hist(samples[, , k], breaks = seq(theta_true[k] - 3, theta_true[k] + 3, 0.1), main = paste("k = ", k),
		col = c_dark, border = c_dark_highlight, xlab = "theta", yaxt = 'n', ylab = "")
	abline(v = theta_true[k], col = "white", lwd = 2)
	abline(v = theta_true[k], col = "black", lwd = 1)
}

## Fit data (hierarchical)
model_sec2.1_hier = cmdstan_model("./normal_hierarchical.stan")

results = model_sec2.1_hier$sample(data = data, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = 1000, iter_sampling = 1000, save_warmup = TRUE,
	max_treedepth = 13, adapt_delta = 0.8)

results$cmdstan_diagnose()

samples = results$draws("tau")
samples[samples > 50] = 50

par(mfrow=c(1, 1))

prior_samples = rgamma(1e4, shape = 0.001, rate = 0.001)
prior_samples = prior_samples[prior_samples < 50]
hist(prior_samples, breaks = seq(0, 50, 0.5),
	main = "", col = c_light, border = c_light_highlight,
	xlab = "tau", yaxt = 'n', ylab = "", ylim = c(0, 300))
hist(samples, breaks = seq(0, 50, 0.5),
	col = c_dark, border = c_dark_highlight, add = TRUE)


samples = results$draws("theta")

par(mfrow=c(3, 3))

for (k in 1:9)
{
	hist(samples[, , k], breaks = seq(theta_true[k] - 3, theta_true[k] + 3, 0.1), main = paste("k = ", k),
		col = c_dark, border = c_dark_highlight, xlab = "theta", yaxt = 'n', ylab = "")
	abline(v = theta_true[k], col = "white", lwd = 2)
	abline(v = theta_true[k], col = "black", lwd = 1)
}

#### Section 2.2: Mixture Population Models
## Common variables
n_chains = 4

## Fit data (discrete mixture model)
model_sec2.2_discrete_mixture = cmdstan_model("./normal_hierarchical_discrete_mixture.stan")

results = model_sec2.2_discrete_mixture$sample(data = data, parallel_chains = n_chains, refresh = 100, chains = n_chains,
	iter_warmup = 1000, iter_sampling = 1000, save_warmup = TRUE, max_treedepth = 13, adapt_delta = 0.8, seed = 4938483)

results$cmdstan_diagnose()

## Divergence ruins the party
theta_array = results$draws("theta")
tau_1_array = results$draws("tau1")
tau_2_array = results$draws("tau2")

divergences = which(results$sampler_diagnostics()[, , "divergent__"] == 1)

par(mfrow=c(3, 3))

for (k in 1:9)
{
	plot(theta_array[, , k], log(tau_1_array), col = c_dark_trans, pch = 16, main = paste("k = ", k),
		xlab = paste0("theta[", k, "]"), xlim = c(theta_true[k] - 3, theta_true[k] + 3),
		ylab = "log(tau1)", ylim = c(-7, 0))

	points(as.vector(theta_array[, , k])[divergences], as.vector(log(tau_1_array))[divergences],
		col = c_green_trans, pch = 16)
}

par(mfrow=c(3, 3))

for (k in 1:9)
{
	plot(theta_array[, , k], log(tau_2_array), col = c_dark_trans, pch = 16, main = paste("k = ", k),
		xlab = paste0("theta[", k, "]"), xlim = c(theta_true[k] - 3, theta_true[k] + 3),
		ylab = "log(tau1)", ylim = c(1.5, 4.5))

	points(as.vector(theta_array[, , k])[divergences], as.vector(log(tau_2_array))[divergences],
		col = c_green_trans, pch = 16)
}

# Why are the funnels manifesting only in the parameters with small true values? Recall that a centred parameterization of an individual
# parameter is optimal when the corresponding individual likelihood function is more informative than the population model, and a
# non-centred parameterization is optimal when the population model is more informative. When the population model is specified by a single
# normal probability density function the dominant contribution is largely determined by the width of the individual likelihood function
# relative to the width of the normal population model, but not its location.

