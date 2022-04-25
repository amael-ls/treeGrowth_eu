functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real totalTreeWeight,
		real averageGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		return (exp(averageGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0>[n_indiv] totalTreeWeight; // Sum of the tree weights for a given plot at a given tieme
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_indiv] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);
	// Normalised and centred totalTreeWeight
}

parameters {
	// Parameters
	real averageGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd(Yobs)^2> sigmaProc;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> sigmaObs; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	// Declaration output variables
	array [n_obs] real newObservations;
	array [n_obs] real latent_dbh_parentsChildren;
	array [n_latentGrowth] real procError_sim;
	array [n_latentGrowth + n_indiv] real yearly_latent_dbh;
	
	{
		// Variables declared in nested blocks are local variables, not generated quantities, and thus won't be printed.
		real current_latent_dbh;
		int growth_counter = 1;
		int dbh_counter = 1;
		real expected_growth;

		for (i in 1:n_indiv) // Loop over all the individuals
		{
			// Fill the parents
			latent_dbh_parentsChildren[parents_index[i]] = latent_dbh_parents[i];
			yearly_latent_dbh[dbh_counter] = latent_dbh_parents[i];

			// Generate the parent observation conditional on the parent state
			newObservations[parents_index[i]] = normal_rng(latent_dbh_parents[i], sigmaObs);

			// Starting point to compute latent dbh child
			current_latent_dbh = latent_dbh_parents[i];

			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1], normalised_totalTreeWeight[i],
					averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);

				// Record difference expected growth versus latent growth
				procError_sim[growth_counter] = expected_growth - latent_growth[growth_counter];

				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter]; // Or should it be += gamma(mean = latent_growth, var = sigmaProc)?
				growth_counter += 1;
				dbh_counter += 1;
				yearly_latent_dbh[dbh_counter] = current_latent_dbh;
			}

			// The last current_latent_dbh corresponds to the child latent dbh
			latent_dbh_parentsChildren[children_index[i]] = current_latent_dbh;
			newObservations[children_index[i]] = normal_rng(latent_dbh_parentsChildren[children_index[i]], sigmaObs);
			dbh_counter += 1;
		}
	}
}