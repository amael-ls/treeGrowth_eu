
/*
	I think it would have been more correct to use the beta binomial distribution for the data. Unfortunately, this distribution is
		unstable when it is reparameterised with mean and variance (rather that its two shape parameters alpha and beta). For instance, with
		the following parameters: N = 10 (the max), alpha = 600, and beta = 400, I get the anlytical mean and variance:
			mean = 6
			variance = 2.421578...

		Now, if I sample 1e3 numbers from the beta bionomial, I get for instance the following mean and variance:
			mean = 5.966
			variance = 2.349193
		Although not too far (we are talking about a difference of around 0.2%), this couple of mean and variance provides the following α
			alpha = -225.3685,
		which is really far from 600; it even makes stan to reject alpha as it must be positive.

		Therefore, I switched to a normal distribution for now (continuous, despite the observations are discrete... Should be changed)
*/

data {
	// Number of data
	int<lower = 1> n_trees; // Number of trees

	// Data
	int<lower = 0> dbh1[n_trees]; // dbh measured by person 1
	int<lower = 0> dbh2[n_trees]; // dbh measured by person 2
	vector<lower = 0>[n_trees] expected_dbh_2016; // Expected dbh
	vector<lower = 0>[n_trees] sd_growth;
}

parameters {
	// Parameter error
	real<lower = 0> error;

	// Latent state
	vector<lower = 0>[n_trees] latent_dbh;
}

model {
	// Defines variables
	target += gamma_lpdf(error | 3.0^2/1, 3.0/1); // Gives an average of 3 and a variance of 1

	for (i in 1:n_trees)
	{
		target += normal_lpdf(latent_dbh[i] | expected_dbh_2016[i], sd_growth[i]);
		target += normal_lpdf(dbh1[i] | latent_dbh[i], error);
		target += normal_lpdf(dbh2[i] | latent_dbh[i], error);
	}
}
