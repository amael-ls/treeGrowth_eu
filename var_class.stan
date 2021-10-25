
data
{
	// Size integer
	int<lower = 1> N; // dimension

	// Vector data
	vector[N] var_class; // Response var
	vector[N] diam_class; // Explanatory variable
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real<lower = 0.0000001> a; // intercept
	real<lower = 0.0000001> b; // decreasing rate
	real<lower = 0.0000001> sigma; // sd
}

model
{
	// Declare variables
	real avg[N];
	for (i in 1:N)
		avg[i] = a/diam_class[i]^b;
		// avg[i] = a * exp(-b * diam_class[i]); // To flat and below variance for dbh > 300

	// Prior
	target += gamma_lpdf(a | 1.0/10000, 1.0/10000);
	target += gamma_lpdf(b | 1.0/10000, 1.0/10000);
	target += gamma_lpdf(sigma | 1.0/10000, 1.0/10000);

	// Model
	target += normal_lpdf(var_class | avg, sigma);
}
