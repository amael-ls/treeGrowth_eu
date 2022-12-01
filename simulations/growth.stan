
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
	real growth(real dbh0, real temperature, real averageGrowth, real dbh_slope, real dbh_slope2,
		real tas_slope, real tas_slope2)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- temperature: The temperature
			
			It takes the following parameters:
				- averageGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- dbh_slope2: The slope for dbh (quadratic term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, quadratic term)
		*/

		return (averageGrowth + dbh_slope*dbh0 + dbh_slope2*dbh0^2 + tas_slope*temperature + tas_slope2*temperature^2);
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

	// Explanatory variables
	real<lower = 0> sd_dbh; // To standardise the initial dbh (sd_dbh is however the sd of all the dbh, not only initial ones)

	array[n_indiv, n_latent_growth] real tas; // Temperature
	real tas_mu; // To centre the temperature
	real<lower = 0> tas_sd; // To standardise the temperature
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
	real averageGrowth;
	real dbh_slope;
	real dbh_slope2;

	real tas_slope;
	real tas_slope2;
	
	// Errors (observation and process)
	// --- Process error, which is the sdlog parameter of a lognormal distrib /!\
	real<lower = 0.5/sd_dbh^2> sigmaProc;
	
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
	target += normal_lpdf(averageGrowth | 0, 20);
	target += normal_lpdf(dbh_slope | 0, 20);
	target += normal_lpdf(dbh_slope2 | 0, 20);

	target += normal_lpdf(tas_slope | 0, 20);
	target += normal_lpdf(tas_slope2 | 0, 20);

	// --- Errors
	target += lognormal_lpdf(sigmaProc | 0.2468601 - log(sd_dbh), 0.16);

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		temporary = latent_dbh_parents[i];
		temporary_tm1 = temporary;
		obs_counter = 1;
		for (j in 1:n_latent_growth) // Loop over growing years
		{
			// Process model
			expected_growth_meanlog = growth(temporary, normalised_tas[i, j],
				averageGrowth, dbh_slope, dbh_slope2, tas_slope, tas_slope2);

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
	
	// Prior on initial hidden state: This is a diffuse initialisation
	target += uniform_lpdf(latent_dbh_parents | 0.1/sd_dbh, 3000/sd_dbh); // Remember that the dbh is in mm and standardised

	// --- Observation model (likelihood)
	for (i in 1:n_indiv)
	{
		for (j in 1:n_obs_growth)
		{
			target += normal_lpdf(normalised_avg_yearly_growth_obs[i, j] | latent_avg_yearly_growth[i, j], 3/sd_dbh);
		}
	}

}
