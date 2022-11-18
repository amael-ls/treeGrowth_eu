
#### Aim of prog: Replicate the small simulation done in Aubry-Kientz et. al (2017) to verify if Eitzel 2013 can be applied for Δt >= 5

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)

#### Simulate data
## Common variables
set.seed(1969-08-18) # Woodstock seed

n = 1000
n_years = 25

init_dbh = rnorm(n = n, mean = 50, sd = 10)
if (any(init_dbh < 0))
	warning("There are some negative starting dbh")

temperature = rnorm(n = n_years, mean = 7, sd = 2.5)

## Parameters (that should be recovered via Stan)
beta0 = 2
beta1 = 0.02
beta2 = 1.5

## Data
treeData = data.table(tree_id = 1:n, dbh1 = init_dbh)
treeData[, paste0("dbh", 2:n_years) := 0]

for (i in 2:n_years)
{
	current_dbh = paste0("dbh", i)
	previous_dbh = paste0("dbh", i - 1)

	treeData[, c(current_dbh) := (1 + beta1)*treeData[, ..previous_dbh] + beta0 + beta2*temperature[i]]
}
