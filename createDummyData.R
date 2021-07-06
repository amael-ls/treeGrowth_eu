
#### Aim of prog: create dummy data for the model designed to parameterise growth
## Description:
#	1. Generate mu_m, the mesured circumference
#	2. Generate mu_r, the real (i.e., true) circumference given there is no error measurement on mu_m
#	3. Generate Y, the real diameter of the tree, that include error measurements
#
## Explanations:
# I generate yearly data, but will parameterise the model with data available only
# every five years. The unknown data between two census must, therefore, been infered.
# 
# Here, we assume that a good year is followed by a bad year with a proba 0.7, and that a bad year is followed
# by a good year with a proba 0.55:
#		P(Good | Bad) = 0.55
#		P(Bad | Bad) = 0.45
#		P(Good | Good) = 0.3
#		P(Bad | Good) = 0.7
#
# With the property:
#	P(Good | Bad) + P(Bad | Bad) = 1
#	P(Good | Good) + P(Bad | Good) = 1
#
# On a good year, the poisson parameter is 8, on a bad year it is 2.
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

library(data.table)

#### Generate data
## Common variables
n = 200 # Number of individuals
if (n %% 2 != 0)
	warning("It would be better to choose an even number")

time = 0:10
id = rep(1:n, each = length(time))

nbMeasurements = length(id)

set.seed(1969-08-18) # Woodstock seed

## Parameters to estimate
# Good times, bad times
goodLambda = 8 # Average growth on good years
badLambda = 2 # Average growth on bad years

# proba transitions
pBtoG = 0.7 # P(Good | Bad) = 0.7
pGtoG = 0.3 # P(Good | Good) = 0.3
pBtoB = 0.55 # P(Bad | Bad) = 0.55
pGtoB = 0.45 # P(Bad | Good) = 0.45

## mu_m
# Store data
data = data.table(tree_id = id, year = rep(time, n), yearQuality = character(nbMeasurements),
	mu_m = integer(nbMeasurements), expectedIncrFromPrevYear = integer(nbMeasurements))

# Initialising, year 0
data[year == 0, mu_m := round(runif(n, 25, 300))]
data[year == 0, yearQuality := sample(x = c("good", "bad"), size = n, replace = TRUE)]
data[year == 0, expectedIncrFromPrevYear := NA]

# Iterating (will be the Hidden Markov process)
for (i in 1:max(time))
{
	# Pick-up lambda. cf comment for proba good/bad years
	previousState = data[year == i - 1, yearQuality]
	currentState = character(n)
	currentState[previousState == "good"] = sample(x = c("good", "bad"), size = sum(previousState == "good"),
		replace = TRUE, prob = c(pGtoG, pGtoB))
	currentState[previousState == "bad"] = sample(x = c("good", "bad"), size = sum(previousState == "bad"),
		replace = TRUE, prob = c(pBtoG, pBtoB))

	data[year == i, yearQuality := currentState]
	data[year == i, expectedIncrFromPrevYear := ifelse(yearQuality == "good", goodLambda, badLambda)]

	data[year == i, mu_m := data[year == i - 1, mu_m] + rpois(n, data[year == i, expectedIncrFromPrevYear])]
}

## Real circumference, the error term is Unif(-1, 1)/2 given it is arounding error from the measuring tape
data[, mu_r := mu_m + runif(nbMeasurements, -1, 1)/2]

## Increments measured and real #* It is assumed that the data table is ordered by year within indiv!
data[, incr_m := mu_m - shift(mu_m, 1), by = tree_id]
data[, incr_r := mu_r - shift(mu_r, 1), by = tree_id]

## Response variable (this will come later, for now let us assume that the resp. var is mu_r)
#! TO DO

#### Check the generated data
## Increment average for good and bad years. Should match goodLambda and badLambda
data[yearQuality == "good", mean(incr_m, na.rm = TRUE)]
data[yearQuality == "bad", mean(incr_m, na.rm = TRUE)]

## Check the probabilities of transitions based on counts
data[, prevYearQuality := shift(yearQuality, 1), by = tree_id]

# From good to good
(estim_pGtoG = data[(yearQuality == "good") & (prevYearQuality == "good"), .N]/data[(yearQuality == "good") & (!is.na(prevYearQuality)), .N])
pGtoG

# From bad to good
(estim_pBtoG = data[(yearQuality == "good") & (prevYearQuality == "bad"), .N]/data[(yearQuality == "good") & (!is.na(prevYearQuality)), .N])
pBtoG

estim_pGtoG + estim_pBtoG

# From good to bad
(estim_pGtoB = data[(yearQuality == "bad") & (prevYearQuality == "good"), .N]/data[(yearQuality == "bad") & (!is.na(prevYearQuality)), .N])
pGtoB

# From bad to bad
(estim_pBtoB = data[(yearQuality == "bad") & (prevYearQuality == "bad"), .N]/data[(yearQuality == "bad") & (!is.na(prevYearQuality)), .N])
pBtoB

estim_pGtoB + estim_pBtoB

#### Available data for the study
## I keep the first and last record for each tree, + a random record in between first and last
kept_indices = sort(
	c(seq(0, nbMeasurements - 1, by = length(time)) + 1, # First indiv measurement
	seq(length(time), nbMeasurements, by = length(time)), # Last indiv measurement
	seq(0, nbMeasurements - 1, by = length(time)) +
		sample(x = 2:(length(time) - 1), size = n, replace = TRUE)))

censusData = data[kept_indices]
