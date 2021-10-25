
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
	real<lower = 0.0001> alpha;
	real<lower = 0.0001> beta;
	real<lower = 0.0001, upper = 3> gamma;
	real<lower = 0.0001> epsilon;
}

model
{
	// Declare variables
	vector[N] expected = 1 + alpha * exp(-beta * dbh);
	vector[N] sigma;
	for (i in 1:N)
	{
		sigma[i] = epsilon/dbh[i]^gamma;
	}

	// Prior
	target += gamma_lpdf(alpha | 0.2/1000, 0.2/1000);
	target += gamma_lpdf(beta | 0.05/1000, 0.05/1000);
	target += gamma_lpdf(gamma | 1.0/10, 1.0/10);
	target += gamma_lpdf(epsilon | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000

	// Model
	for (i in 1:N)
	{
		target += gamma_lpdf(Y[i] | expected[i]^2/sigma[i], expected[i]/sigma[i]);
	}
}
