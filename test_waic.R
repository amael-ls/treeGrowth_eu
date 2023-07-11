
library(data.table)
library(cmdstanr)
library(loo)
options(mc.cores = 4)

n = 5e3

intercept = -5.35
slope = 1.325
sigma = 2.9

x = runif(n, min = -1, max = 12)
y = rnorm(n = n, mean = intercept + slope*x, sd = sigma)

stanData = list(
	n = n,
	x = x,
	y = y
)

model = cmdstan_model("test_waic.stan")

results = model$sample(data = stanData, parallel_chains = 4, refresh = 500, chains = 4)

mean(results$draws("intercept"))
mean(results$draws("slope"))
mean(results$draws("sigma"))

loglik = results$draws("log_lik")
dim(loglik)

n_chains = results$num_chains()
n_iter = results$metadata()$iter_sampling

waic = function(loglik_draws, n_chains, n_iter)
{
	S = n_chains*n_iter
	
	mean_draws = apply(X = exp(loglik_draws), MARGIN = 3, FUN = mean) # The exp is necessary because stan takes the log of proba distrib fct
	lpd_hat = sum(log(mean_draws))

	n_indiv_rw = length(mean_draws)
	p_hat = numeric(n_indiv_rw)

	for (measure in 1:n_indiv_rw)
	{
		temporary = (loglik_draws[, , measure] - 1/S*sum(loglik_draws[, , measure]))^2
		p_hat[measure] = 1/(S - 1) * sum(temporary)

	}
	print(paste("min = ", min(p_hat)))
	print(paste("max = ", max(p_hat)))

	return (c(lpd_waic = lpd_hat, p_waic = sum(p_hat), waic = -2*(lpd_hat - sum(p_hat))))
}

waic(loglik_draws = loglik, n_chains = n_chains, n_iter = n_iter)

loo::waic(loglik)
