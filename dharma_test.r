
#### Ain of prog: Check if I am able to use DHARMa correctly
## Plan of script:
#	1. Create data
#	2. Fit the data with stan
#	3. Check the residuals

#### Load libraries and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(cmdstanr)
library(DHARMa)

options(max.print = 500)

#### Create data
## Common variables
n_data = 2500
min_x = -2
max_x = 10

## Parameters
a = 2.5
b = 0.78
sigma = 6

## Data
x = runif(n = n_data, min = min_x, max = max_x)
y = rnorm(n = n_data, mean = a + b*x, sd = sigma)

#### Fit data
model = cmdstan_model("dharma_test.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = n_data,
	x = x,
	y = y)

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

#### Plot
avg = function(x, a = a_ptw, b = b_ptw)
	return (a + b*x);

a_ptw = mean(results$draws(variables = c("a"), format = "draws_df")$a)
b_ptw = mean(results$draws(variables = c("b"), format = "draws_df")$b)
sigma_ptw = mean(results$draws(variables = c("sigma"), format = "draws_df")$sigma)

plot(x, y, pch = 19, col = "#323232AA")
curve(avg(x), col = "#CD212A", from = min_x, to = max_x, lwd = 2, add = TRUE)

#### Residuals
n_rep = 250

dt = data.table(rep_x = rep(x, each = n_rep), sampled = numeric(n_rep * length(x)))

for (i in 1:length(x))
{
	avg_i = avg(x[i]);
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rnorm(n = n_rep, mean = avg_i, sd = sigma_ptw)]
}

sims = matrix(data = dt[, sampled], nrow = n_rep, ncol = length(x)) # each column is for one data point (the way I organised the dt)
sims = t(sims) # Transpose the matrix for dharma

forDharma = createDHARMa(simulatedResponse = sims, observedResponse = y)
plot(forDharma)
