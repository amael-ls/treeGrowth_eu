
/*
	I decided to reparameterise the gamma distribution in function of mean and variance rather than shape and rate. Remember that we have:
		shape = mean^2/var
		rate = mean/var
*/

data {
	// Number of data
	int<lower = 1> n_trees; // Number of trees

	// Data
	vector<lower = 0>[n_trees] growth;
}

parameters {
	real<lower = 0> mu; // mean
	real<lower = 0> sigma; // variance (and not sd!)
}

model {
	// Defines variables
	target += gamma_lpdf(mu | 3.0^2/5.0, 3.0/5.0);
	target += gamma_lpdf(sigma | 5.0^2/100.0, 5.0/100.0);
	target += gamma_lpdf(growth | mu^2/sigma, mu/sigma);
}
