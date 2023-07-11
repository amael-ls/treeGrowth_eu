
data {
	int<lower = 1> n;

	vector[n] x;
	vector[n] y;
}

parameters {
	real intercept;
	real slope;
	real<lower = 0> sigma;
}

model {
	// Priors
	target += normal_lpdf(intercept | 0, 10);
	target += normal_lpdf(slope | 0, 10);
	target += gamma_lpdf(sigma | 1.0/100.0, 1.0/100.0);

	// Likelihood
	for (i in 1:n)
		target += normal_lpdf(y[i] | intercept + slope*x[i], sigma);
}

generated quantities {
	vector[n] log_lik;
	for (i in 1:n)
		log_lik[i] = normal_lpdf(y[i] | intercept + slope*x[i], sigma);
}

