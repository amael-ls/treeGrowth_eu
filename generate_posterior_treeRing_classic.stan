
// This stand-alone script generates annual growth with the SSM approach for a given time-series of environmental conditions and an initial diameter

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
	int<lower = 1> n_climate_new; // Dimension of environmental vectors
	int<lower = 1> time_series_length; // Number of new latent growth = sum of the length of each time series (for tree rings data)
	array[n_indiv_new] int<lower = 1, upper = time_series_length - 1> start_ind; // Indexing for environment
	array[n_indiv_new] int<lower = 1> nbGrowingYears; // Individual time series length (for growth, dbh would need a +1 to each compound)

	array [n_indiv_new] real dbh0; // Initial dbh (which is not dbh_init from growth.stan)

	array [n_climate_new] real pr_new; // precipitation, already standardised
	array [n_climate_new] real tas_new; // temperature, already standardised
	array [n_climate_new] real ph_new; // pH, already standardised
	array [n_climate_new] real basal_new; // basal area, already standardised
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
	array [time_series_length] real simulatedGrowth; // Simulated latent annual growth
	array [time_series_length] real simulatedGrowth_avg; // Averaged latent annual growth (i.e., average of a lognormal)
	
	array [time_series_length + n_indiv_new] real simulatedLatentDBH; // Simulated latent annual dbh
	array [time_series_length + n_indiv_new] real simulatedLatentDBH_avg; // Averaged latent annual dbh (i.e., from the averaged growth)
	array [n_indiv_new] real simulatedObservedDBH; // Simulated last observation per individual
	array [n_indiv_new] real simulatedObservedDBH_avg; // Simulated last observation per individual
	
		{
		// Common variables, any variable defined here will not be part of the output
		int count_dbh = 1;
		int count_growth = 1;
		int j;
		array [4] real env;
		real meanlog = 0;
		real meanlog_avg = 0;
		real current_dbh = 0;
		real current_dbh_avg = 0;
		real sigmaObs = 3/sd_dbh; // Error on dbh, do not use 2*sigmaObs/delta_years which is for growth
		
		// Generate growth for different climate and dbh combinations
		for (i in 1:n_indiv_new)
		{
			current_dbh = dbh0[i]/sd_dbh;
			current_dbh_avg = dbh0[i]/sd_dbh;

			simulatedLatentDBH[count_dbh] = current_dbh;
			simulatedLatentDBH_avg[count_dbh] = current_dbh_avg;

			// for (j in start_ind[i]:end_ind[i])
			for (k in 1:nbGrowingYears[i])
			{
				j = start_ind[i] + k - 1; // The -1 is because k should be from 0 to nbYears - 1
				// Extract the environment in the following order ---> precip, temperature, pH and competition
				env = {pr_new[j], tas_new[j], ph_new[j], basal_new[j]};

				// Compute growth
				meanlog = growth(current_dbh, env, averageGrowth, dbh_slope, dbh_slope2,
					pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
				meanlog_avg = growth(current_dbh_avg, env, averageGrowth, dbh_slope, dbh_slope2,
					pr_slope, pr_slope2 , tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
			
				simulatedGrowth[count_growth] = lognormal_rng(meanlog, sigmaProc);
				current_dbh += simulatedGrowth[count_growth]; // This is latent dbh!

				simulatedGrowth_avg[count_growth] = exp(meanlog_avg + sigmaProc^2/2);
				current_dbh_avg += simulatedGrowth_avg[count_growth]; // This is averaged latent dbh!

				count_dbh += 1;
				count_growth += 1;
				simulatedLatentDBH[count_dbh] = current_dbh;
				simulatedLatentDBH_avg[count_dbh] = current_dbh_avg;
			}
			simulatedObservedDBH[i] = normal_rng(current_dbh, sigmaObs);
			simulatedObservedDBH_avg[i] = normal_rng(current_dbh_avg, sigmaObs);
			count_dbh += 1; // The growth counter stops one year earlier! Indeed n measurements implies only n - 1 growing years!
		}
	}
}
