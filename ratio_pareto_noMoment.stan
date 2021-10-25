
data
{
	// Size integer
	int<lower = 1> N; // Size of response var

	// Vector data
	vector[N] Y; // Response var
	vector[N] dbh; // Explanatory variable
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real<lower = 0.00000001> alpha_0;
	real<lower = 0.00000001> alpha_1;
}

model
{
	// Declare variables
	vector[N] alpha;
	for (i in 1:N)
	{
		alpha[i] = alpha_0*dbh[i]^alpha_1; // alpha MUST be an increasing function of dbh, to decrease the average and variance
	}

	// Prior
	target += gamma_lpdf(alpha_0 | 1.0/100000, 1.0/100000);
	target += gamma_lpdf(alpha_1 | 1.0/100000, 1.0/100000);

	// Model
	for (i in 1:N)
	{
		target += pareto_lpdf(Y[i] | 1, alpha[i]);
	}
}
