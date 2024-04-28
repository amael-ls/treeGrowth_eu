
// This stand-alone script generates posterior median for a fixed dbh, stand basal area and pH, and for many climatic values (from a raster)

functions {
	// Function for growth. This returns the expected log(growth) in mm, for 1 year, i.e. the parameter meanlog of lognormal distribution
	real growth(real dbh0, array[] real x_r, real averageGrowth, real dbh_slope, real dbh_slope2,
		real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real ph_slope, real ph_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- x_r: array of real data in the following order ---> precip, temperature, pH and competition
			
			It takes the following parameters:
				- averageGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- dbh_slope2: The slope for dbh (quadratic term)
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (quadratic term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, quadratic term)
				- ph_slope: The slope for pH
				- ph_slope2: The slope for pH (quadratic term)
				- competition_slope: The slope for competition
		*/

		real precip = x_r[1];
		real temperature = x_r[2];
		real ph = x_r[3];
		real standBasalArea = x_r[4];

		return (averageGrowth + dbh_slope*dbh0 + dbh_slope2*dbh0^2 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + ph_slope*ph + ph_slope2*ph^2 + competition_slope*standBasalArea);
	}
}

data {
	// Old data required for the parameters
	int<lower = 1> n_indiv; // Total number of individuals (all NFIs together)
	int<lower = 1> n_obs; // Total number of tree observations (all NFIs together)
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_growth; // Number of measured growth = n_obs - n_indiv
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data
	real<lower = 0> sd_dbh; // To standardise the lower and upper bounds of dbh

	// New data
	int<lower = 0> N_draws; // Number of draws to compute the probability
	real dbh0; // Initial dbh (which is not dbh_init from growth.stan)
	array [4] real env_current; // Contains in this order: pr, tas, ph, basal area; pH and basal area are fixed, while pr and tas come from current raster
	array [4] real env_future; // Contains in this order: pr, tas, ph, basal area; pH and basal area are fixed, while pr and tas come from raster future climate
}

parameters {
	// Parameters for growth function
	real averageGrowth;
	real dbh_slope;
	real dbh_slope2;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real ph_slope;
	real ph_slope2;

	real competition_slope;
	
	// Errors (observation and process)
	// --- Process error, which is the sdlog parameter of a lognormal distrib /!\
	real<lower = 0.5/sd_dbh^2, upper = 10> sigmaProc;

	// --- Extreme error, by default at least twice the min observation error. Rüger 2011 found it is around 8 times larger
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> proba; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
	vector<lower = 0>[n_growth] latent_avg_annual_growth; // Real (and unobserved) averaged annual growth
}

generated quantities {
	real median_current = 0;
	real median_future = 0;
	real<lower = 0, upper = 1> probaGrowth_inf_median = 0;
	{
		// Common variables, any variable defined here will not be part of the output
		real sigmaObs = 3.0/sd_dbh; // Fixed value from growth.stan
		real meanlog_future = 0;
		
		array [N_draws] int is_simulated_growth_lower;

		// Compute meanlog
		median_current = exp(growth(dbh0, env_current, averageGrowth, dbh_slope, dbh_slope2,
			pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope));

		meanlog_future = growth(dbh0, env_future, averageGrowth, dbh_slope, dbh_slope2,
			pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);

		median_future = exp(meanlog_future);

		for (i in 1:N_draws)
			is_simulated_growth_lower[i] = lognormal_rng(meanlog_future, sigmaProc) < median_current;

		probaGrowth_inf_median = mean(is_simulated_growth_lower); // Compute the proba
	}
}
