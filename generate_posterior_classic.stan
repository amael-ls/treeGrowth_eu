
#### This stand-alone script generates annual growth with the classic approach for a given time-series of environmental conditions and an initial diameter

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

	int<lower = 1> n_climate_new; // Dimension of the new climate vector
	int<lower = n_climate_new - 1, upper = n_climate_new - 1> n_years; // Number of growing years, it can only be n_climate_new - 1

	real dbh0; // Initial dbh (which is not dbh_init from growth.stan)

	array [4*n_climate_new] real x_r; // Contains in this order: pr, ts, ph, basal area, each of size n_climate_new, and already standardised
	array [4] real x_r_avg; // Contains in this order: pr, ts, ph, basal area averaged over n_climate_new years
}

transformed data {
	real normalised_dbh0 = dbh0/sd_dbh;
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
	real simulatedGrowth_avg_clim_avg;
	real final_dbh;
	array [n_years] real simulatedGrowth_avg;
	array [n_climate_new] real current_dbh;
	current_dbh[1] = normalised_dbh0;
	{
		// Common variables, any variable defined here will not be part of the output
		array [4] int index;
		real meanlog = 0;

		// Generate growth for different climate and dbh combinations
		for (i in 1:n_years)
		{
			// Extract the i^th value of precipitation, temperature, pH, and basal area, respectively
			index = {i, i + n_climate_new, i + 2*n_climate_new, i + 3*n_climate_new};
			meanlog = growth(current_dbh[i], x_r[index], averageGrowth, dbh_slope, dbh_slope2,
				pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
			
			simulatedGrowth_avg[i] = exp(meanlog + sigmaProc^2/2.0);
			current_dbh[i + 1] = current_dbh[i] + simulatedGrowth_avg[i];
		}
		meanlog = growth(normalised_dbh0, x_r_avg, averageGrowth, dbh_slope, dbh_slope2,
				pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
		simulatedGrowth_avg_clim_avg = n_years*exp(meanlog + sigmaProc^2/2.0);
		final_dbh = normalised_dbh0 + simulatedGrowth_avg_clim_avg;
	}
}
