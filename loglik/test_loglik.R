
#### Aim of prog: Fit a simple model with stan, and compute manually its log-likelihood

#### Simulate data
n = 1e4
mu = 0.4*x
sigma = 3.7

set.seed(1969-08-18) # Woodstock seed

x = runif(n = n, min = 0, max = 25)
y = rnorm(n = n, mean = mu, sd = sigma)

#### Compute minus log-likelihood
## Priors
prior_mu = function(mu)
	return (- dnorm(x = mu, mean = 0, sd = 100, log = TRUE));

prior_sigma = function(sigma)
	return (- dgamma(x = sigma, shape = 5^2/100, rate = 5/100, log = TRUE)); # i.e., mean = 5, var = 100

## Log-likelihood
loglik = function(data, mu, sigma)
	return (- sum(dnorm(x = data, mean = mu, sd = sigma, log = TRUE)));

## Log probability lp__, WATCH OUT: the minus signs are in the priors and loglik functions already
lp = function(data, mu, sigma)
	return (loglik(data, mu, sigma) + prior_mu(mu) + prior_sigma(sigma));

#### Fit data with stan model

