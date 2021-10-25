
data {
	int<lower = 1> nb_data;

	real init_state;

	real observations[nb_data];
}

parameters {
	real states[nb_data];

	real<lower = 0> sigmaObservation;
	real<lower = 0> sigmaProcess;
}

model {
	// Priors
	target += gamma_lpdf(sigmaObservation | 1.0/10, 1.0/10);
	target += gamma_lpdf(sigmaProcess | 5^2/10.0, 5.0/10);

	// States
	target += normal_lpdf(states[1] | init_state, sigmaProcess);

	for (i in 2:nb_data)
		target += normal_lpdf(states[i] | states[i - 1], sigmaProcess);

	// Likelihood
	target += normal_lpdf(observations | states, sigmaObservation);
}
