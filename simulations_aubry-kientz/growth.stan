
/*
	Comments:
		- Yobs is an vector, containing all dbh observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
		- The observation error for the French data is different from the other countries because there have been corrections in the data.
			Therefore, there is no extrem error due to, for instance, typos.
		- Note that the dbh is standardised but not centred. The reason was that in a previous version of the model, it would have
			implied a division by 0. Since then, I do not have this problem anymore, but I kept the dbh non-centred.
	
	Reminder (method of first and second moments, i.e., mean and variance or standard deviation):
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
		
		- The lognormal distribution of stan uses, in this order, meanlog (mu) and sdlog (sigma). It can however be reparametrised
			using mean and sd: mu = log(mean^2/sqrt(sd^2 + mean^2)), sigma = sqrt(log(sd^2/mean^2 + 1))

	Like C++, BUGS, and R, Stan uses 0 to encode false, and 1 to encode true. https://mc-stan.org/docs/functions-reference/logical-functions.html
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real temperature, real beta0, real beta1, real beta2,
		real beta3, real beta4)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- temperature: The temperature
			
			It takes the following parameters:
				- beta0: The basic growth, when all the explanatory variables are set to zero
				- beta1: The slope for dbh
				- beta2: The slope for dbh (quadratic term)
				- beta3: The slope for temperature (tas)
				- beta4: The slope for temperature (tas, quadratic term)
		*/

		return (beta0 + beta1*dbh0 + beta2*dbh0^2 + beta3*temperature + beta4*temperature^2);
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Total number of individuals
	int<lower = 1> n_obs_growth; // Total number of growth observations
	int<lower = 1> n_latent_growth; // Number of growing years (for latent state)
	int<lower = 1> delta_t; // Number of growing years (for latent state)

	// Observations
	array[n_indiv, n_obs_growth] real avg_yearly_growth_obs;
	array[n_indiv] real <lower = 0> dbh_init; // Initial dbh data

	// Explanatory variables
	real<lower = 0> sd_dbh; // To standardise the initial dbh (sd_dbh is however the sd of all the dbh, not only initial ones)

	array[n_indiv, n_latent_growth] real tas; // Temperature
	real tas_mu; // To centre the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	// Provided parameters
	// real beta0;
	// real beta1;
	// real beta2;
	// real beta3;
	// real beta4;

	// Observation error
	real <lower = 0> sigmaObs;
}

transformed data {
	// Normalised but NOT centred averaged yearly growth
	array[n_indiv, n_obs_growth] real normalised_avg_yearly_growth_obs = avg_yearly_growth_obs;
	for (i in 1:n_indiv)
		for (j in 1:n_obs_growth)
	 		normalised_avg_yearly_growth_obs[i, j] /= sd_dbh;

	// Normalised and centred temperatures
	array[n_indiv, n_latent_growth] real normalised_tas;
	for (i in 1:n_indiv)
		for (j in 1:n_latent_growth)
	 		normalised_tas[i, j] = (tas[i, j] - tas_mu)/tas_sd;
}

parameters {
	// Parameters for growth function
	real beta0;
	real beta1;
	real beta2;

	real beta3;
	real beta4;
	
	// Errors (observation and process)
	// --- Process error, which is the sdlog parameter of a lognormal distrib /!\
	real<lower = 0, upper = 2> sigmaProc;
	
	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
	array[n_indiv, n_latent_growth] real<lower = 0> latent_growth; // Real (and unobserved) yearly growth
}

model {
	// Declare variables
	int obs_counter; // Counter in the 'observation space' for chidlren
	real expected_growth_meanlog;
	real temporary;
	real temporary_tm1; // temporary at time t - 1 (useful for trees measured more than twice)
	array[n_indiv, n_obs_growth] real latent_avg_yearly_growth;
	
	// Priors
	// --- Growth parameters
	target += normal_lpdf(beta0 | 0, 20);
	target += normal_lpdf(beta1 | 0, 20);
	target += normal_lpdf(beta2 | 0, 20);

	target += normal_lpdf(beta3 | 0, 20);
	target += normal_lpdf(beta4 | 0, 20);

	// --- Errors
	// target += lognormal_lpdf(sigmaProc | 0.2468601 - log(sd_dbh), 0.16);
	target += uniform_lpdf(sigmaProc | 0, 2);

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		// Initial values
		temporary = latent_dbh_parents[i];
		temporary_tm1 = temporary;
		obs_counter = 1;
		
		// Prior on initial hidden state: remember that latent_dbh_parents is in mm and standardised.
		target += gamma_lpdf(latent_dbh_parents[i] | dbh_init[i]^2/5.0, sd_dbh*dbh_init[i]/5.0);

		// Markov process
		for (j in 1:n_latent_growth) // Loop over growing years
		{
			// Process model
			expected_growth_meanlog = growth(temporary, normalised_tas[i, j],
				beta0, beta1, beta2, beta3, beta4);

			target += lognormal_lpdf(latent_growth[i, j] | expected_growth_meanlog, sigmaProc);

			// Dbh at time t + 1
			temporary += latent_growth[i, j];

			// Only the relevant (i.e., children) diameters at breast height are recorded
			// print("j = ", j);
			if (j % delta_t == 0)
			{
				// print("(from if) j = ", j);
				latent_avg_yearly_growth[i, obs_counter] = (temporary - temporary_tm1)/delta_t;
				obs_counter += 1;
				temporary_tm1 = temporary;
			}
		}
	}

	// --- Observation model (likelihood)
	for (i in 1:n_indiv)
	{
		for (j in 1:n_obs_growth)
		{
			target += normal_lpdf(normalised_avg_yearly_growth_obs[i, j] | latent_avg_yearly_growth[i, j], 2*sigmaObs/delta_t);
		}
	}

}
