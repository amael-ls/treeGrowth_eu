
// To encode our sparsity assumption we'll use a beta(3,1) prior density function for γ that allocates 0.875 of its total probability
//	above γ = 0.5.
// This sentence can be understood this way: integrate(dbeta, lower = 0.5, upper = 1, shape1 = 3, shape2 = 1)
//	0.875 with absolute error < 9.7e-15

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
	real<lower = 0, upper = 1> gamma;
	real<lower = 0> tau1;
	real<lower = 0> tau2;
}

model {
	for (k in 1:K)
		target += log_mix(gamma, normal_lpdf(theta[k] | 0, tau1), normal_lpdf(theta[k] | 0, tau2));
	
	gamma ~ beta(3, 1);
	
	tau1 ~ normal(0, 0.1);
	tau2 ~ normal(0, 10);

	// Observational model
	y ~ normal(theta[context_idx], sigma);
}
