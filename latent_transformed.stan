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
	real<lower = 0.0001> latentMean;
	real<lower = 0.0001> latentVar;
	real<lower = 0.0001> measurementError;

	vector<lower = 0>[N] trueDbh;
}

transformed parameters
{
	real latentShape = latentMean^2/latentVar;
	real latentRate = latentMean/latentVar;
}

model
{
	// Prior
	target += gamma_lpdf(latentMean | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000
	target += gamma_lpdf(latentVar | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000
	target += gamma_lpdf(measurementError | 0.0001, 0.0001); // Gives an average of 1 and variance of 10000

	// Absolute determinant of the jacobian transformation
	target += 2*log(latentMean) - log(latentVar);

	// Increment target log probability density with gamma_lupdf(y | alpha, beta)
	target += gamma_lpdf(trueDbh | latentShape, latentRate);
	target += normal_lpdf(Y | trueDbh, measurementError);
}
