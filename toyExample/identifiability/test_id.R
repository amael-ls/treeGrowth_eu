
library(bayesplot)
library(cmdstanr)

#### Generate data
## Common variables
sigmaObservation = 2
sigmaProcess = 8

nb_data = 2000

states = numeric(nb_data)
Y = numeric(nb_data)

set.seed(1969-08-18) # Woodstock seed

## Process
# Initial state
states[1] = runif(1, 10, 15)

# States
for (i in 2:nb_data)
	states[i] = rnorm(1, states[i - 1], sigmaProcess)

## Observations
Y = rnorm(nb_data, states, sigmaObservation) # Equivalent to Y = Id * states + rnorm(nb_data, 0, sigmaObservation)

#### Estimate parameters
## Data
init_state = rnorm(1, states[1], 1)
stanData = list(
	# Number of data
	nb_data = nb_data,
	init_state = init_state,
	observations = Y
)

## Compile model
model = cmdstan_model("./test_id.stan")

## Run model
# Common variables
n_chains = 4
maxIter = 2000

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 500, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2)

results$cmdstan_diagnose()

plot_title = ggplot2::ggtitle("Posterior distribution of sigmaObservation", "with medians and 80% intervals")
mcmc_areas(results$draws("sigmaObservation"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigmaObservation, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for sigmaObservation")
mcmc_trace(results$draws("sigmaObservation")) + plot_title

plot_title = ggplot2::ggtitle("Posterior distribution of sigmaProcess", "with medians and 80% intervals")
mcmc_areas(results$draws("sigmaProcess"), prob = 0.8) + plot_title +
	ggplot2::geom_vline(xintercept = sigmaProcess, color = "#FFA500")

plot_title = ggplot2::ggtitle("Traces for sigmaProcess")
mcmc_trace(results$draws("sigmaProcess")) + plot_title
