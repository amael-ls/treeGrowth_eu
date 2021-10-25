
data
{
	// Size integer
	int<lower = 1> N; // Size of response var

	// Vector data
	vector[N] x; // Explanatory variable
	vector[N] y; // Response variable
}

parameters
{
	real a;
	real b;
	real<lower = 0.000001> sigma;
}

model
{
	// Prior
	target += normal_lpdf(a | 0, 10000);
	target += normal_lpdf(b | 0, 10000);
	target += gamma_lpdf(sigma | 1.0/100000, 1.0/100000);

	target += normal_lpdf(y | a + b*x, sigma);
}
