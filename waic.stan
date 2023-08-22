
/*
	This stan script is too compute the WAIC of both models (SSM and Classic approaches). The WAIC computation is the same between the two
	because the growth processus is the same. Only the estimated parameters theta_G differ!

	Note that I ended all the variables related to the ring-width time series by rw (for ring-width)
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real ph, real standBasalArea, real averageGrowth, real dbh_slope, real dbh_slope2,
		real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real ph_slope, real ph_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temperature: The temperature
				- ph: The pH (CaCl2)
				- standBasalArea: Basal area of the stand (i.e., includes all the trees, regardless of the species)
			
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

		return (averageGrowth + dbh_slope*dbh0 + dbh_slope2*dbh0^2 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + ph_slope*ph + ph_slope2*ph^2 + competition_slope*standBasalArea);
	}
}

data {
	// Dimensions
	int<lower = 1> n_data_rw; // Number of ring width data

	// Explanatory variables
	// --- Diameters
	vector[n_data_rw] dbh_rw; //! Already standardised
	real<lower = 0> sd_dbh;
	
	// --- Climate
	vector[n_data_rw] precip_rw; //! Already standardised
	vector[n_data_rw] tas_rw; //! Already standardised
	
	// --- Soil
	vector[n_data_rw] ph_rw; //! Already standardised
	
	// --- Basal area
	vector[n_data_rw] standBasalArea_rw; //! Already standardised
	
	// Observed data
	vector[n_data_rw] ring_width; // Correspond to annual growth,  //! Already standardised
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

	// No need to put the other parameters (latent states and observation error)
}

generated quantities {
	// Declared quantities here are part of the output
	vector[n_data_rw] log_lik;

	{
		// Declared quantities in this block will not be returned
		int climate_counter = 1;
		real meanlog;

		// Compute log-likelihood
		for (i in 1:n_data_rw)
		{
			meanlog = growth(dbh_rw[i], precip_rw[i],
				tas_rw[i], ph_rw[i], standBasalArea_rw[i],
				averageGrowth, dbh_slope, dbh_slope2, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
			
			log_lik[i] = lognormal_lpdf(ring_width[i] | meanlog, sigmaProc); // note that _lpdf already gives the log proba density function
		}
	}
}

