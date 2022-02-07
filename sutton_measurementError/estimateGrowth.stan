
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
	// real<lower = 0> mu; // mean
	real a; // mean = a dbh^2 + b^dbh + c
	real b; // mean = a dbh^2 + b^dbh + c
	real c; // mean = a dbh^2 + b^dbh + c
	real<lower = 0> sigma; // variance (and not sd!)
}

model {
	// Defines variables
	vector[n_trees] mu;

	mu = exp(a*dbh_standardised^2 + b*dbh_standardised + c);

	// Priors
	target += normal_lpdf(b | 0, 10); // variance of 10 is huge, given that dbh is std and that mu is on exponential!
	target += normal_lpdf(a | 0, 10); // variance of 10 is huge, given that dbh is std and that mu is on exponential!
	target += normal_lpdf(c | 0, 10); // variance of 10 is huge, given that dbh is std and that mu is on exponential!
	target += gamma_lpdf(sigma | 5.0^2/100.0, 5.0/100.0);

	// Model
	// target += gamma_lpdf(mu | 3.0^2/5.0, 3.0/5.0);
	target += gamma_lpdf(growth | mu^2/sigma, mu/sigma);
}
