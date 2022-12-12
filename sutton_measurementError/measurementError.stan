
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

functions {
	// Expected growth for a given dbh
	real mu_G(real dbh0, real g, real h, real i, real scaling){
		return (g*exp(-exp(h - i*dbh0/scaling)));
	}

	// Expected variance of growth for a given dbh
	real var_G(real dbh0, real a, real b, real c, real scaling){
		return (exp(a*(dbh0/scaling)^2 + b*dbh0/scaling + c));
	}
}

data {
	// Number of data
	int<lower = 1> n_trees; // Number of trees (individuals)
	int<lower = 1> n_states; // Number of states

	// indices
	array [n_trees] int<lower = 1> index_parents; // index of the measurement corresponding to the first campaign (2011-2013)
	array [n_trees] int<lower = 1> index_children; // index of the measurement corresponding to the campaign in 2016

	// Parameters growth
	real a; // For var_G
	real b; // For var_G
	real c; // For var_G

	real<lower = 0> g; // For mu_G
	real h; // For mu_G
	real i; // For mu_G

	real<lower = 0> scaling; // Scaling dbh that was USED to parameterise growth

	// Data
	array [n_trees] int<lower = 0> dbh0; // dbh measured in the first campaign (2011-2013) neither by person 1 nor 2
	array [n_trees] int<lower = 0> dbh1; // dbh measured by person 1 in 2016
	array [n_trees] int<lower = 0> dbh2; // dbh measured by person 2 in 2016
	array [n_trees] int<lower = 0> years; // Number of years spent between the first campaign and 2016
}

parameters {
	// Parameter error
	real<lower = 0> error;

	// Latent state
	vector<lower = 0>[n_states] latent_dbh; // Latent dbh from first campaign to 2016 for all the individuals
}

model {
	// Defines variables
	real expected_dbh;
	real var_growth;

	// Computation with EM algorithm (compute condition Expectations, and Maximise likelihood with respect to params)
	for (indiv in 1:n_trees)
	{
		// Expected dbh a year after the first campaign
		expected_dbh = latent_dbh[index_parents[indiv]] + mu_G(latent_dbh[index_parents[indiv]], g, h, i, scaling);
		var_growth = var_G(latent_dbh[index_parents[indiv]], a, b, c, scaling);

		for (ind in (index_parents[indiv] + 1):(index_children[indiv])) // The first being the parent is accounted for later.
		{
			target += gamma_lpdf(latent_dbh[ind] | (expected_dbh)^2/var_growth, expected_dbh/var_growth);
			expected_dbh = latent_dbh[ind] + mu_G(latent_dbh[ind], g, h, i, scaling);
			var_growth = var_G(latent_dbh[ind], a, b, c, scaling);
		}
	// End EM algo
		// Likelihood parents
		target += normal_lpdf(dbh0[indiv] | latent_dbh[index_parents[indiv]], error);

		// Likelihood children (measured by person 1 and 2)
		target += normal_lpdf(dbh1[indiv] | latent_dbh[index_children[indiv]], error);
		target += normal_lpdf(dbh2[indiv] | latent_dbh[index_children[indiv]], error);
	}

	// Priors
	target += gamma_lpdf(error | 3.0^2/1, 3.0/1); // Gives an average of 3 and a variance of 1
	target += lognormal_lpdf(latent_dbh[index_parents] | 5.246511, 0.5269077); // cf measurementError.R for these values
}
