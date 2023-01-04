
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
	vector<lower = 0>[n_trees] dbh0;
}

transformed data {
	vector[n_trees] dbh_standardised = dbh0/sd(dbh0);
}

parameters {
	real a;
	real b;
	real c;

	real<lower = 0> g;
	real h;
	real i;
}

model {
	// Defines variables
	vector[n_trees] mu;
	vector[n_trees] sigma;
	mu = g*exp(-exp(h - i*dbh_standardised));
	sigma = exp(a*dbh_standardised^2 + b*dbh_standardised + c);

	// Priors
	target += normal_lpdf(a | 0, 10);
	target += normal_lpdf(b | 0, 10);
	target += normal_lpdf(c | 0, 10);

	target += gamma_lpdf(g | 3.0^2/10, 3.0/10);
	target += normal_lpdf(h | 0, 100);
	target += normal_lpdf(i | 0, 100);
	
	// Model
	for (ind in 1:n_trees)
		target += gamma_lpdf(growth[ind] | mu[ind]^2/sigma[ind], mu[ind]/sigma[ind]);
}
