
#### Aim of prog: create dummy data for a hierarchical model

#### Clear memory and load libraries
rm(list = ls())
graphics.off()

library(data.table)
library(viridis)

#### Create data
## Common variables
set.seed(1969-08-18) # Woodstock seed

n_group = 4

min_per_group = 25
max_per_group = 40

data_per_group = sample(x = min_per_group:max_per_group, size = n_group, replace = TRUE)
n_data = sum(data_per_group)

colours = viridis(n_group)

## Parameters
beta0 = 5.5
sigma_beta = 2
beta0_group = rnorm(n = n_group, mean = beta0, sd = sigma_beta)

beta1 = 2.3

sigma_res = sqrt(2)

realParameters = c(beta0, sigma_beta, beta0_group, beta1, sigma_res)
names(realParameters) = c("beta0", "sigma_beta", paste0("beta0_", 1:n_group), "beta1", "sigma_res")

## Explanatory variable
x = runif(n = n_data, min = -5, max = 10)
dt = data.table(y = numeric(n_data), x = x, group = rep(1:n_group, data_per_group))

dt[, colour := colours[group], by = group]
dt[, .N, by = group] # This should be equal to data_per_group

count = 0
for (i in 1:n_group)
{
	ind = (1 + count):(data_per_group[i] + count)
	dt[ind, y := rnorm(n = data_per_group[i], mean = beta0_group[i] + beta1*x, sd = sigma_res)]
	count = count + data_per_group[i]
}

plot(dt[, x], dt[, y], pch = 19, col = dt[, colour])

#### Save data
saveRDS(dt, "./data.rds")
saveRDS(realParameters, "realParameters.rds")
