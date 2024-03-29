
/*
	Comments:
		- avg_annual_growth_obs is a vector, containing all growth observations (averaged per year, i.e., (dbh1 - dbh0)/(t1 - t0))
		- I call "parent" the first measured dbh of each individual.
		- Note that the dbh is standardised but not centred. The reason was that in a previous version of the model, it would have
			implied a division by 0. Since then, I do not have this problem anymore, but I kept the dbh non-centred.
		- The error structure is as follow:
			1. sigmaProc (process error) is the unexplained variance for growth
			2. sigmaObs (observation error) is the routine error done when measuring dbh
			3. etaObs (observation error) is the extreme error done when measuring dbh occurring with probability proba
		- Note that the routine observation error, sigmaObs, is a fixed parameter taken from the literature. Indeed it was impossible to
			estimate this parameter together with etaObs, proba, and sigmaProc.
		- Never forget that it does not make sense to directly compare the parameters estimated with the State-Space Model vs the classic
			approach done here! Although the parameters have the same names between those from the SSM and those from the classic approach,
			they do not have the same meaning. Example: pr_slope is the responce of annual growth to annual precipitation in SSM, while it
			is the response of averaged growth to averaged precipitations over delta_t years in the classic approach.
	
	Reminder (method of first and second moments, i.e., mean and variance or standard deviation, see Appendix C 'Rescaling'):
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var

		- When the gamma distribution is scaled by a constant C, it is equivalent to multiply the rate by C
		
		- The lognormal distribution of stan uses, in this order, meanlog (mu) and sdlog (sigma). It can however be reparametrised
			using mean and var: mu = log(mean^2/sqrt(var + mean^2)), sigma = sqrt(log(var/mean^2 + 1))

		- When the lognoraml distribution is scaled by a constant C, it is equivalent to substract log(C) to the meanlog parameter

	Like C++, BUGS, and R, Stan uses 0 to encode false, and 1 to encode true. https://mc-stan.org/docs/functions-reference/logical-functions.html
*/
functions {
	// Function for growth. This returns the expected log(growth) in mm, for 1 year, i.e., meanlog parameter of lognormal distribution
	real avg_growth(real dbh0, real precip, real temperature, real ph, real standBasalArea, real averageGrowth, real dbh_slope,
		real dbh_slope2, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real ph_slope, real ph_slope2,
		real competition_slope)
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
	// Number of data
	int<lower = 1> n_indiv; // Total number of individuals (all NFIs together)
	int<lower = 1> n_climate; // Dimension of the climate vector (all NFIs together)
	int<lower = 1, upper = n_indiv> n_plots; // Number of plots (all NFIs together)
	int<lower = 1> n_obs; // Total number of tree observations (all NFIs together)
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_growth; // Number of measured growth = n_obs - n_indiv
	array [n_indiv] int<lower = 1, upper = n_obs> nbIntervalGrowth; // Number of growth intervals for each individual
	array [n_growth] int<lower = 1> deltaYear; // Number of years between two measurements of an individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

	// Indices
	array [n_indiv] int<lower = 1, upper = n_climate> climate_index; // Index of the climate associated to each growth interval
	
	array [n_inventories] int<lower = 1, upper = n_growth - 1> start_nfi_avg_growth; // Starting point of each NFI for averaged obs growth
	array [n_inventories] int<lower = 2, upper = n_growth> end_nfi_avg_growth; // Ending point of each NFI for averaged obs growth
	
	array [n_indiv] int<lower = 1, upper = n_plots> plot_index; // Indicates to which plot individuals belong to

	// Observations
	vector[n_growth] avg_annual_growth_obs;
	vector<lower = 0>[n_indiv] dbh_init; // Initial dbh data

	// Explanatory variables
	real<lower = 0> sd_dbh; // To standardise the initial dbh (sd_dbh is however the sd of all the dbh, not only initial ones)

	vector<lower = 0>[n_climate] precip; // Precipitations
	real<lower = 0> pr_mu; // To centre the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To centre the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0, upper = 14>[n_plots] ph; // pH of the soil measured with CaCl2
	real<lower = 0, upper = 14> ph_mu; // To centre the pH
	real<lower = 0> ph_sd; // To standardise the pH

	vector<lower = 0>[n_climate] standBasalArea; // Sum of the tree basal area for a given plot at a given time (interpolation for NA data)
	real ba_mu; // To centre the basal area
	real<lower = 0> ba_sd; // To standardise the basal area
}

transformed data {
	vector[n_growth] normalised_avg_annual_growth_obs = avg_annual_growth_obs/sd_dbh; // Normalised but NOT centred averaged annual growth
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	vector[n_plots] normalised_ph = (ph - ph_mu)/ph_sd; // Normalised and centred pH
	vector[n_climate] normalised_standBasalArea = (standBasalArea - ba_mu)/ba_sd; // Normalised and centred BA
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

model {
	// Declare variables
	int growth_counter = 1; // Counter in the latent space
	real expected_avg_growth_meanlog;
	real current_latent_dbh; // current_latent_dbh at time t

	// --- Routine observation error, Assessing Precision in Conventional Field Measurements of Individual Tree Attributes (Luoma 2017)
	real sigmaObs = 3.0/sd_dbh; // Std. Dev. of a normal distrib /!\
	
	// Priors
	// --- Growth parameters
	target += normal_lpdf(averageGrowth | 0, 20);
	target += normal_lpdf(dbh_slope | 0, 20);
	target += normal_lpdf(dbh_slope2 | 0, 20);
	
	target += normal_lpdf(pr_slope | 0, 20);
	target += normal_lpdf(pr_slope2 | 0, 20);

	target += normal_lpdf(tas_slope | 0, 20);
	target += normal_lpdf(tas_slope2 | 0, 20);

	target += normal_lpdf(ph_slope | 0, 20);
	target += normal_lpdf(ph_slope2 | 0, 20);

	target += normal_lpdf(competition_slope | 0, 20);

	// --- Errors
	/*
		A note on the folowing priors:
		--> Let m be the mean of the gamma distribution, and v its variance (see comment at the beginning of this file to get shape and
			rate from mean and variance).

		For instance, for the routine measure error, which is the standard deviation of the observation around the latent state, we have:
		m = sqrt(3)
		v = 0.025

		Then, because I work on standardised dbh, I need to divide m by sd(dbh), and v by sd(dbh)^2. Therefore, we get the following:
		shape = (m/sd(dbh))^2 / (v/sd(dbh)^2) = m^2/v
		rate =  (m/sd(dbh))   / (v/sd(dbh)^2) = sd(dbh) * m/v

		--> Let m be the mean of the lognormal distribution, and s its standard deviation (see comment at the beginning of this file to get
		meanlog and sdlog from mean and sd).
		For the process error, which is the sd of growth, we have:
		meanlog = log( (m/sd(dbh))^2/sqrt((s/sd(dbh))^2 + (m/sd(dbh))^2) ) = log(m^2/sqrt(s^2 + m^2)) - log(sd(dbh)) = meanlog - log(sd(dbh))
		sdlog = sqrt(log( (s/sd(dbh))^2 / (m/sd(dbh))^2 + 1)) = sqrt(log(s^2/m^2 + 1))

		The values are taken from Rüger et al (2011), Growth Strategies of Tropical Tree Species: Disentangling Light and Size Effects
		for etaObs and proba priors.

		The values are taken from Luoma et al (2017), Assessing Precision in Conventional Field Measurements of Individual Tree Attributes
		for sigmaObs prior.
	*/
	target += uniform_lpdf(sigmaProc | 0.5/sd_dbh^2, 10); // I suppose that the process error is between 0 and 10 (quite large, mean and var involve Exp[sdlog^2]!)
	target += gamma_lpdf(etaObs | 30^2/45.0, sd_dbh*30/45.0); // <=> extreme measurement error (sd) = 30 mm ± 6.7 mm
	target += beta_lpdf(proba | 48.67, 1714.84); // This corresponds to a 2.76 % chance extrem error, ± 0.39 % (Rüger et. al 2011)

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		// Prior on initial hidden state: remember that latent_dbh_parents is in mm and standardised.
		target += gamma_lpdf(latent_dbh_parents[i] | dbh_init[i]^2/5.0, sd_dbh*dbh_init[i]/5.0);

		current_latent_dbh = latent_dbh_parents[i];
		for (j in 1:nbIntervalGrowth[i]) // Loop over number of intervals measured for individual i
		{
			// Process model
			expected_avg_growth_meanlog = avg_growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
			normalised_tas[climate_index[i] + j - 1], normalised_ph[plot_index[i]], normalised_standBasalArea[climate_index[i] + j - 1],
			averageGrowth, dbh_slope, dbh_slope2, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);

			target += lognormal_lpdf(latent_avg_annual_growth[growth_counter] | expected_avg_growth_meanlog, sigmaProc);

			// New dbh at time t + Δt: dbh1 = dbh0 + (t1 - t0) * avg_annual_growth
			current_latent_dbh += deltaYear[growth_counter]*latent_avg_annual_growth[growth_counter];

			growth_counter += 1;
		}
	}
	
	// --- Observation model
	for (k in 1:n_inventories)
	{
		/*
			Compare true (i.e., hidden or latent) latent averaged annual growth with observed averaged annual growth
			Do not try to vectorise here! https://mc-stan.org/docs/2_29/stan-users-guide/vectorizing-mixtures.html

			Note that here, I use 2 * error/Δt. This is because the annual growth error is twice the error on the dbh,
			divided by the number of years between the two measurement (i.e., averageing the growth error). See Rüger et al 2011.
		*/
		for (i in start_nfi_avg_growth[k]:end_nfi_avg_growth[k])
			target += log_mix(proba[k],
				normal_lpdf(normalised_avg_annual_growth_obs[i] | latent_avg_annual_growth[i], 2*etaObs[k]/deltaYear[i]),
				normal_lpdf(normalised_avg_annual_growth_obs[i] | latent_avg_annual_growth[i], 2*sigmaObs/deltaYear[i]));
	}
}
