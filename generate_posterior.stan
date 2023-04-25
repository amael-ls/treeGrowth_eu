
// This stand-alone script generates annual growth series with the SSM approach for given time-series of environmental conditions and initial diameters

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
	int<lower = 1> n_indiv_new; // Number of new individuals (which corresponds to the number of initial diameters provided)
	int<lower = 1> n_latent_dbh_new; // Number of years (for dbh)
	int<lower = n_latent_dbh_new - n_indiv_new, upper = n_latent_dbh_new - n_indiv_new> n_growing_years_new; // Number of growing years (new)
	array [n_indiv_new] int<lower = 1, upper = n_latent_dbh_new> nbYearsGrowth_new; // Number of years of growth for new individuals

	array [n_indiv_new] real dbh0; // Initial dbh (which is not dbh_init from growth.stan)

	array [4*n_latent_dbh_new] real x_r; // Contains in this order: pr, tas, ph, basal area, each of size n_latent_dbh_new, and already standardised
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
	real<lower = 0.5/sd_dbh^2> sigmaProc;

	// --- Extreme error, by default at least twice the min observation error. Rüger 2011 found it is around 8 times larger
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> proba; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
	vector<lower = 0>[n_growth] latent_growth; // Real (and unobserved) yearly growth
}

generated quantities {
	array [n_growing_years_new] real simulatedGrowth;
	array [n_growing_years_new] real simulatedGrowth_avg;
	array [n_indiv_new] real simulatedObservedGrowth;
	
	array [n_latent_dbh_new] real current_dbh;
	array [n_latent_dbh_new] real current_dbh_obs;
	array [n_latent_dbh_new] real current_dbh_avg;
	
	{
		// Common variables, any variable defined here will not be part of the output
		array [4] int index;
		real meanlog = 0;
		real meanlog_avg = 0;
		real sigmaObs = 3.0/sd_dbh; // Fixed value from growth.stan

		int dbh_count = 1;
		int growth_count = 1;
		int origin_time_series_i;
		
		// Generate growth for different climate and dbh combinations
		for (i in 1:n_indiv_new)
		{
			current_dbh[dbh_count] = dbh0[i]/sd_dbh;
			current_dbh_obs[dbh_count] = normal_rng(current_dbh[dbh_count], sigmaObs);
			current_dbh_avg[dbh_count] = dbh0[i]/sd_dbh;

			origin_time_series_i = dbh_count;

			for (j in 1:nbYearsGrowth_new[i])
			{
				// Extract the i^th value of precipitation, temperature, pH, and basal area, respectively
				index = {dbh_count, dbh_count + n_latent_dbh_new, dbh_count + 2*n_latent_dbh_new, dbh_count + 3*n_latent_dbh_new};

				// Compute growth
				meanlog = growth(current_dbh[dbh_count], x_r[index], averageGrowth, dbh_slope, dbh_slope2,
					pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
				meanlog_avg = growth(current_dbh_avg[dbh_count], x_r[index], averageGrowth, dbh_slope, dbh_slope2,
					pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
			
				simulatedGrowth[growth_count] = lognormal_rng(meanlog, sigmaProc);
				simulatedGrowth_avg[growth_count] = exp(meanlog_avg + sigmaProc^2/2);
				
				current_dbh[dbh_count + 1] = current_dbh[dbh_count] + simulatedGrowth[growth_count]; // This is latent dbh!
				current_dbh_obs[dbh_count + 1] = normal_rng(current_dbh[dbh_count], sigmaObs); // This is what could have been potentially measured
				current_dbh_avg[dbh_count + 1] = current_dbh_avg[dbh_count] + simulatedGrowth_avg[growth_count]; // This is averaged latent dbh!

				dbh_count += 1;
				growth_count += 1;
			}

			// Generate potential observation, dbh error converted to growth error by dividing by the time interval
			simulatedObservedGrowth[i] = normal_rng(current_dbh[dbh_count] - current_dbh[origin_time_series_i], 2*sigmaObs/nbYearsGrowth_new[i]);
			
			dbh_count += 1;
		}
	}
}
