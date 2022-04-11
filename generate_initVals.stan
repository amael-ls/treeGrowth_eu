
/*
	Comments:
		- Yobs is an vector, containing all dbh observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
		- The observation error for the French data is different from the other countries because there have been corrections in the data.
			Therefore, there is no extrem error due to, for instance, typos.
		- Note that the dbh is standardised but not centred. The reason was that in a previous version of the model, it would have
			implied a division by 0. Since then, I do not have this problem anymore, but I kept the dbh non-centred.
	
	Reminder:
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real totalTreeWeight,
		real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temperature: The temperature
				- totalTreeWeight: The sum of all the tree weights
			
			It takes the following parameters:
				- potentialGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, squared term)
				- competition_slope: The slope for competition
		*/

		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_indiv_new; // New number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = 1> n_latentGrowth_new; // Dimension of the new state space vector for latent growth
	int<lower = n_obs - n_indiv_new, upper = n_obs - n_indiv_new> n_children; // Number of children trees observations = n_obs - n_indiv_new
	array [n_indiv_new] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv_new] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv_new] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	/*
		The tree weight correspond to the weight of the tree per hectar. It accounts for:
			1. the area of the plot (disc with a radius of 6, 9 or 15 metres depending on the size class)
			2. the proximity to the boundary (rectification of the probability of drawing trees at the edge)
	*/
	vector<lower = 0>[n_indiv_new] totalTreeWeight; // Sum of the tree weights for a given plot at a given tieme
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	vector[n_indiv_new] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);
	// Normalised and centred totalTreeWeight
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd_dbh^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd_dbh> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0.1/sd_dbh, upper = 2000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	// Declare variables
	vector [n_latentGrowth_new] latent_growth_init;
	array [n_indiv_new] real latent_dbh_parents_init = normal_rng(normalised_Yobs[parents_index], 0.5/sd(Yobs)); // Start close to the measure

	{
		real latent_dbh_children;
		real expected_growth;
		int growth_counter = 1;

		// Model
		for (i in 1:n_indiv_new) // Loop over all the individuals
		{
			latent_dbh_children = latent_dbh_parents_init[i];
			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(latent_dbh_children, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1], normalised_totalTreeWeight[i],
					potentialGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);
				latent_growth_init[growth_counter] = gamma_rng(expected_growth^2/processError, expected_growth/processError);

				// Dbh at time t + 1, only the last (i.e., child) dbh is recorded
				latent_dbh_children += latent_growth_init[growth_counter];
				growth_counter += 1;
			}
		}
	}
}
