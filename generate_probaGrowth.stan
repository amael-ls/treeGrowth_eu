
/*
	This script computes the probability that annual yearly growth is above a threshold, given the estimated parameters
	with the data used to estimate them, and new data (i.e., new predictors).
*/

functions {
	// Function computing the expected annual diameter-growth in mm for a given dbh, environment, and set of parameters
	real growth(real dbh0, real precip, real temperature, real ph, real standBasalArea, real averageGrowth, real dbh_slope,
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
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, squared term)
				- ph_slope: The slope for pH
				- ph_slope2: The slope for pH (squared term)
				- competition_slope: The slope for competition
		*/

		return (exp(averageGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + ph_slope*ph + ph_slope2*ph^2 + competition_slope*standBasalArea));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1, upper = n_indiv> n_plots; // Number of plots (all NFIs together)
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual
	array [n_children] int<lower = 1> deltaYear; // Number of years between two measurements of an individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

	int<lower = 1> n_dbh; // Dimension of the dbh vector (new data)
	int<lower = 1> n_threshold; // Dimension of the threshold vector (new data, lower bound of integral, i.e., variable of cdf function)

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

	vector<lower = 0>[n_dbh] dbh0; // DBH used as a predictor (new data)
	vector<lower = 0>[n_threshold] threshold; // DBH used as a predictor (new data)

	array [4] real environment; // Contains in this order: precip, temp, pH, basal area, already standardised
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_dbh] normalised_DBH = dbh0/sd_dbh; // Normalised but NOT centred dbh
	vector[n_threshold] normalised_threshold = threshold/sd_dbh; // Normalised threshold
	vector[n_children] normalised_avg_yearly_growth_obs = avg_yearly_growth_obs/sd_dbh; // Normalised but NOT centred averaged yearly growth
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centered temperatures
	vector[n_plots] normalised_ph = (ph - ph_mu)/ph_sd; // Normalised and centred pH
	vector[n_climate] normalised_standBasalArea = (standBasalArea - ba_mu)/ba_sd; // Normalised and centred BA
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
	array[n_dbh*n_threshold] real probaGrowth_beyondThreshold;
	{
		real expected_growth;
		real growth_mean_logNormal;
		real growth_sd_logNormal;
		int count = 1;

		// for (i in 1:n_dbh)
		for (i in 1:n_threshold)
		{
			for (j in 1:n_dbh)
			{
				expected_growth = growth(normalised_DBH[j],
					environment[1], environment[2], environment[3], environment[4],
					averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
				
				growth_mean_logNormal = log(expected_growth^2/sqrt(sigmaProc^2 + expected_growth^2));
				growth_sd_logNormal = sqrt(log(sigmaProc^2/expected_growth^2 + 1));

				// Compute the probability that growth > threshold for a given dbh, environment, and set of parameters.
				probaGrowth_beyondThreshold[count] = 1 - lognormal_cdf(normalised_threshold[i] | growth_mean_logNormal, growth_sd_logNormal);

				count += 1;
			}
		}
	}	
}
