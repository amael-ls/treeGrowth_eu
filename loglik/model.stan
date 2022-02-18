
data {
	int<lower = 1> N; // Number of observations

	vector[N] x;

	vector[N] obs;
}

parameters {
	real slope;
	real<lower = 0> sigma;
}

model {
	// Declare variable
	vector[N] mu = slope * x;

	// Priors
	target += normal_lpdf(slope | 0, 100);
	target += gamma_lpdf(sigma | 5.0^2/100.0, 5.0/100.0);

	// Likelihood
	target += normal_lpdf(obs | mu, sigma);
}