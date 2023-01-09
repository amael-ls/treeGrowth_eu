
/*
	This model is taken from Climate Impacts on Tree Growth in the Sierra Nevada (Aubry-Kientz, 2017). It is based on simulated data
	and tests the effect of different measuring time interval on recovering parameters.
*/

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_years; // Number of years
	int<lower = 1> n_obs_years; // Number of years of diameter observations (=> n_obs_years - 1 growth observations)

	// Observations
	vector<lower = 0> [n_indiv] dbh_init; // Initial diameters
	array[n_indiv, n_obs_years - 1] real <lower = 0> observed_averaged_annual_growth; // Observed averaged annual growth
	int<lower = 1> freq_obs; // Frequency of observation

	// Predictors
	vector[n_years] temperature;
}

parameters {
	// Regression coefficients
	real beta0;
	real beta1;
	real beta2;

	// Latent yearly growth
	array[n_indiv, n_years - 1] real <lower = 0, upper = 100> latent_yearly_growth;

	// Process error
	real<lower = 0> sigmaProc;
}

model {
	// Declaration variables
	real current_dbh;
	real previous_recorded_dbh;
	real current_expected_growth;
	real latent_averaged_annual_growth;

	int counter;

	// Priors
	target += normal_lpdf(beta0 | 0, 20);
	target += normal_lpdf(beta1 | 0, 20);
	target += normal_lpdf(beta2 | 0, 20);

	target += gamma_lpdf(sigmaProc | 1.0/1000.0, 1.0/1000.0);

	// Latent model and likelihood
	for (i in 1:n_indiv)
	{
		previous_recorded_dbh = dbh_init[i];
		current_dbh = dbh_init[i];

		counter = 1;

		for (j in 1:(n_years - 1)) // j is the year growth
		{
			current_expected_growth = beta0 + beta1*current_dbh + beta2*temperature[j];
			target += normal_lpdf(latent_yearly_growth[i, j] | current_expected_growth, sigmaProc); // Process
			
			current_dbh += latent_yearly_growth[i, j];

			if (j % freq_obs == 0) // Likelihood of data
			{
				latent_averaged_annual_growth = (current_dbh - previous_recorded_dbh)/freq_obs;
				target += normal_lpdf(observed_averaged_annual_growth[i, counter] | latent_averaged_annual_growth, sigmaProc);
				previous_recorded_dbh = current_dbh;
				counter = counter + 1;
			}
		}
	}
}



