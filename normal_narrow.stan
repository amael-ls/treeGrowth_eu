
data {
	// Number of data
	int<lower=0> K; // Number of contexts
	int<lower=0> N; // Number of observations

	// Index
	array [N] int<lower = 1, upper = K> context_idx; // Context assignments

	// Std. dev measurement variablility
	real<lower = 0> sigma; // Measurement variability

	// Observations
	vector[N] y;
}

parameters {
	// Horseshoe parameters
	vector[K] theta;
}

model {
	// Horseshoe prior model
	theta ~ normal(0, 0.1);

	// Observational model
	y ~ normal(theta[context_idx], sigma);
}
