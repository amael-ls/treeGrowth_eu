/*
	Comments:
	Gamma distribution uses shape and rate in stan
*/

data
{
	// Size integer
	int<lower = 1> N; // size of response var

	// Vector data
	vector[N] Y; // response var, integers
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real<lower = 1, upper = 3> latentShape;
	real<lower = 0.01, upper = 3> latentRate;
	real<lower = 0.0001> measurementError;

	vector<lower = 0>[N] trueDbh;
}

model
{
	// Prior
	target += uniform_lpdf(latentShape | 1, 5);
	// target += beta_lpdf(latentRate | 0.75, 0.75); // Gives an average of 0.5 and variance of 0.1
	target += uniform_lpdf(latentRate | 0.001, 5);
	// target += gamma_lpdf(measurementError | 4, 2); // Gives an average of 2 and variance of 1
	target += gamma_lpdf(measurementError | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000

	// Increment target log probability density with gamma_lupdf(y | alpha, beta)
	target += gamma_lpdf(trueDbh | latentShape, latentRate);
	target += normal_lpdf(Y | trueDbh, measurementError);
}

generated quantities
{
	real latentMean = latentShape/latentRate;
	real latentVar = latentMean/latentRate;
}
