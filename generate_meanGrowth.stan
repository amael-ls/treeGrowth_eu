
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
	int<lower = 1, upper = n_indiv> n_plots; // Number of plots (all NFIs together)
	int<lower = 1> n_climate_new; // Dimension of the new climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual
	array [n_children] int<lower = 1> deltaYear; // Number of years between two measurements of an individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the 'observation space'
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the 'observation space'
	array [n_children] int<lower = 2> latent_children_index; // Index of children in the 'latent space'
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	array [n_inventories] int<lower = 1, upper = n_indiv - 1> start_nfi_parents; // Starting point of each NFI for parents
	array [n_inventories] int<lower = 2, upper = n_indiv> end_nfi_parents; // Ending point of each NFI for parents
	array [n_inventories] int<lower = 1, upper = n_children - 1> start_nfi_children; // Starting point of each NFI for children
	array [n_inventories] int<lower = 2, upper = n_children> end_nfi_children; // Ending point of each NFI for children
	array [n_inventories] int<lower = 1, upper = n_children - 1> start_nfi_avg_growth; // Starting point of each NFI for averaged obs growth
	array [n_inventories] int<lower = 2, upper = n_children> end_nfi_avg_growth; // Ending point of each NFI for averaged obs growth
	
	array [n_indiv] int<lower = 1, upper = n_plots> plot_index; // Indicates to which plot individuals belong to

	// Observations
	vector<lower = 0>[n_obs] Yobs;
	vector[n_children] avg_yearly_growth_obs;

	real lower_bound;
	real upper_bound;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0, upper = 14>[n_plots] ph; // pH of the soil measured with CaCl2
	real<lower = 0, upper = 14> ph_mu; // To standardise the pH
	real<lower = 0> ph_sd; // To standardise the pH

	vector<lower = 0>[n_climate] standBasalArea; // Sum of the tree basal area for a given plot at a given time (interpolation for NA data)
	real ba_mu; // To standardise the basal area
	real<lower = 0> ba_sd; // To standardise the basal area

	array [3*n_climate_new] real x_r; // Contains in this order: pr, ts, basal area, each of size n_climate_new, and already standardised
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_children] normalised_avg_yearly_growth_obs = avg_yearly_growth_obs/sd_dbh; // Normalised but NOT centred averaged yearly growth
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_plots] normalised_ph = (ph - ph_mu)/ph_sd; // Normalised and centred pH
	vector[n_climate] normalised_standBasalArea = (standBasalArea - ba_mu)/ba_sd; // Normalised and centred BA

	real normalised_lower_bound = lower_bound/sd_dbh;
	real normalised_upper_bound = upper_bound/sd_dbh;

	array [0] int x_i; // Unused, but necessary, int array to compute integral
}

parameters {
	// Parameters for growth function
	real averageGrowth;
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real ph_slope;
	real ph_slope2;

	real competition_slope;
	
	// Errors (observation and process)
	// --- Process error, which is the sd of a lognormal distrib /!\
	real<lower = 0.5/sd_dbh^2> sigmaProc;

	array [n_inventories] real<lower = 0.1/sqrt(12)*25.4/sd_dbh> sigmaObs; // Std. Dev. of a normal distrib /!\
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	array [n_inventories] real<lower = 0, upper = 1> proba; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	array [n_climate_new] real average_growth_sim;
	
	{
		array [3] int index;

		for (i in 1:n_climate_new)
		{
			index = {i, i + n_climate_new, i + 2*n_climate_new}; // To get the i^th value of precip, temp, and basal area, respectively
			average_growth_sim[i] = integrate_1d(growth_integrand, normalised_lower_bound, normalised_upper_bound,
				{averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope},
				x_r[index], x_i)/(normalised_upper_bound - normalised_lower_bound);
		}
	}
}
