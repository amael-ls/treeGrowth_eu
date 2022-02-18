
#### Aim of prog: Fit a simple model with stan, and compute manually its log-likelihood

#### Packages
rm(list = ls())
graphics.off()

library(data.table)
library(cmdstanr)

options(max.print = 500)

#### Simulate data
## Common variables
set.seed(1969-08-18) # Woodstock seed
n = 1e4
x = runif(n = n, min = 0, max = 25)

## Parameters to estimate
slope = 0.4
sigma = 3.7

## Data
mu = slope*x
y = rnorm(n = n, mean = mu, sd = sigma)

data = data.table(x = x, y = y)

#### Compute minus log-likelihood
## Priors
prior_slope = function(slope)
	return (dnorm(x = slope, mean = 0, sd = 100, log = TRUE));

prior_sigma = function(sigma)
	return (dgamma(x = sigma, shape = 5^2/100, rate = 5/100, log = TRUE)); # i.e., mean = 5, var = 100

## Log-likelihood
loglik = function(data, slope, sigma)
	return (sum(dnorm(x = data[, y], mean = slope*data[, x], sd = sigma, log = TRUE)));

## Log probability lp__, WATCH OUT: the minus signs are in the priors and loglik functions already
lp = function(data, slope, sigma)
	return (loglik(data, slope, sigma) + prior_slope(slope) + prior_sigma(sigma) + log(sigma));

# Remark: Why is there a + log(sigma) in the lp function? This is to compensate the transformation (change of variable) that occurs in the
#	back stage of Stan. Indeed, sigma is a constrained variable (sigma > 0), so that a change of variable is used to have it from -Inf to +Inf
#	This implies to multiply by the determinant of the jacobian (or add the log since it is loglik!).
#	See https://vasishth.github.io/bayescogsci/book/sec-firststan.html for more informations

## Compute the average of lp
compute_lp_mean = function(results, data)
{
	n_chains = results$num_chains()
	n_sample = results$metadata()[["iter_sampling"]]

	slope_draws = results$draws("slope")
	sigma_draws = results$draws("sigma")

	lp_vec = numeric(n_chains * n_sample)
	count = 0

	for (col in 1:n_chains)
	{
		for (row in 1:n_sample)
		{
			count = count + 1
			lp_vec[count] = lp(data = data, slope = as.numeric(slope_draws[row, col, 1]), sigma = as.numeric(sigma_draws[row, col, 1]))
		}
	}
	return (mean(lp_vec));
}

#### Fit data with stan model
## Create stan data
stanData = list(
	# Number of data
	N = length(y),

	# Explanatory variable
	x = x,

	# Observations
	obs = y
)

## Compile model
model = cmdstan_model("./model.stan")

## Fit data
n_chains = 4
results = model$sample(data = stanData, parallel_chains = n_chains, chains = n_chains,
	iter_warmup = 1000, iter_sampling = 1000, seed = 1969-08-18)

results$print()

#### Log-likelihood
## Check likelihood from function versus from model
mean(results$draws("lp__"))
compute_lp_mean(results, data)

## Plot log-likelihood and log probability
vec_slope = seq(slope - 0.02, slope + 0.02, length.out = 50)
vec_sigma = seq(sigma - 0.1, sigma + 0.1, length.out = 30)

ll = matrix(data = 0, nrow = length(vec_sigma), ncol = length(vec_slope))
r = 1
for (sig in vec_sigma)
{
	c = 1
	for (sl in vec_slope)
	{
		ll[r, c] = lp(data = data, slope = sl, sigma = sig) # (r - 1)*ncol(ll) + c
		c = c + 1
	}
	r = r + 1
}

x_sol = mean(results$draws("slope"))
y_sol = mean(results$draws("sigma"))
x_real = slope
y_real = sigma

filled.contour(x = vec_slope, y = vec_sigma, z = t(ll), xlab = "Slope", ylab = "Sigma",
	plot.axes = {
		axis(1); axis(2);
		points(x_sol, y_sol, col = "#1D1F54", pch = 19, cex = 2); text(x_sol, y_sol, pos = 3, labels = "Estimated")
		points(x_real, y_real, col = "#1D7FDF", pch = 19, cex = 2); text(x_real, y_real, pos = 1, labels = "Real solution")
	})
