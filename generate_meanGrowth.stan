
functions {
	// Function to integrate (growth function). This returns the expected growth in mm, for 1 year.
	real growth_integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i)
	{
		// theta is the vector of parameters
		real averageGrowth = theta[1];
		real dbh_slope = theta[2];

		real pr_slope = theta[3];
		real pr_slope2 = theta[4];
		real tas_slope = theta[5];
		real tas_slope2 = theta[6];

		real competition_slope = theta[7];

		// x_r is an array of real data (here, precip, competition, etc...)
		real precip = x_r[1];
		real temperature = x_r[2];
		real totalTreeWeight = x_r[3];

		return (exp(averageGrowth + dbh_slope*x + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_climate_new; // Dimension of the new climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual
	array [n_children] int<lower = 1> deltaYear; // Number of years between two measurements of an individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

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

	vector<lower = 0>[n_indiv] totalTreeWeight; // Sum of the tree weights for a given plot at a given time

	array [3*n_climate_new] real x_r; // Contains in this order: pr, ts, totalTreeWeight, each of size n_climate_new, and already standardised
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_indiv] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);

	array [0] int x_i; // Unused, but necessary, int array to compute integral
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

	real<lower = 0.5/sd(Yobs)^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	array [n_climate_new] real average_growth_sim;
	
	{
		real lower_bound = 50/sd(Yobs);
		real upper_bound = 700/sd(Yobs);
		array [3] int index;

		for (i in 1:n_climate_new)
		{
			index = {i, i + n_climate_new, i + 2*n_climate_new}; // To get the i^th value of precip, temp, and tree weight, respectively
			average_growth_sim[i] = integrate_1d(growth_integrand, lower_bound, upper_bound,
				{averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope},
				x_r[index], x_i)/(upper_bound - lower_bound);
		}
	}
}
