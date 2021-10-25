
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
	real<lower = 1, upper = 2> a; // average
	real<lower = 0.0001> b; // average
	real<lower = 0.0001> gamma;
	real<lower = 0.0001> epsilon;
}

model
{
	// Declare variables
	real var_i;
	real avg_i;
	vector[N] lambda;
	vector[N] alpha;
	for (i in 1:N)
	{
		var_i = epsilon/dbh[i]^gamma; // According to var_class.stan, this represent the variance the best
		// avg_i = a * exp(-b*dbh[i]);
		avg_i = a/dbh[i]^b;
		lambda[i] = (1 - avg_i)*(1 - 2*avg_i + avg_i^2 + var_i)/(1 - 2*avg_i + avg_i^2 - var_i); // scale
		alpha[i] = 2*var_i/(-1 + 2*avg_i - avg_i^2 + var_i); // shape
	}

	// Prior
	target += uniform_lpdf(a | 1, 2);
	// target += gamma_lpdf(a | 1.0/1000, 1.0/1000);
	target += gamma_lpdf(b | 1.0/1000, 1.0/1000);
	target += gamma_lpdf(gamma | 1.0/1000, 1.0/1000);
	target += gamma_lpdf(epsilon | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000

	// Model
	for (i in 1:N)
	{
		target += pareto_type_2_lpdf(Y[i] | 1, lambda[i], alpha[i]);
	}
}
